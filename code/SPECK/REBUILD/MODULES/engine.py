"""The device engine: grid + drives + tiers + cavity + the ledger, one BDF.

Composition (all data flows one way, no module-level state):

  params/species/network  ->  CompiledNetwork (chemistry.py)
  loki_tier / kinetic_tier -> electron-impact rates
  beam.py                 ->  deposition profile + conserving yields
  laser.py                ->  cavity photon equation, ledgered loss split
  ledger.py               ->  the S3 closure invariant

STATE: Y = [cells x (species densities, U_e, T_g)] ++ [n_ph] ++ [4 ledger
channels]. One scipy BDF integrates everything, with a block-diagonal
Jacobian sparsity pattern (each cell couples only to itself; the global
photon density couples both ways; ledger rows carry no claimed entries —
they are pure integrals and converge chord-Newton style). That single-solver
choice is what makes the ledger's closure exact to integrator tolerance: the
channel integrals ride the same step controller as the physics.

GRID SEMANTICS: the device grid is COARSE (mm–cm cells; chemistry, not
Debye physics). Cells are independent 0-D reactors coupled through (a) the
shared cavity photon field along z and (b) the 3-D beam deposition profile.
No transport between cells in v1 — at 100–200 ns and cm cells, neutral/ion
transport is negligible; stated here, gated nowhere (it is a scope, not a
claim). The regime line every run prints says exactly this.
"""

from __future__ import annotations

import math
import time as _time

import numpy as np

from . import beam as B
from .chemistry import CompiledNetwork
from .circuit import IE, II1, II2, N_CIRCUIT_STATES, CircuitSpec
from .laser import Cavity, gain_m1, loss_rates, photon_rhs
from .ledger import CHANNELS, N_CHANNELS, closure_report
from .network import build_network, mode_optics
from .params import ParameterSet, engine_parameters
from .species import (ARF_MODE_SPECIES, KRF_MODE_SPECIES, StateLayout)
from .units import (ATM_PA, EV, K_B, M_E, gas_density_m3, photon_energy_J)
from .units import E_CHARGE as E_CHARGE_SIGNED

_UNSET = object()

DEFAULT_MIXTURES = {
    "krf": {"Ar": 0.949, "Kr": 0.05, "F2": 0.001},     # DE-AR Ar-dominant
    "arf_chen": {"Ne": 0.964, "Ar": 0.035, "F2": 0.001, "Xe": 0.0001},
}


def _mw_envelope(t_s: float, t0_ns: float, rise_ns: float,
                 flat_ns: float) -> float:
    """Trapezoidal 0->1 window (the BeamSpec pulse shape, time-shifted)."""
    t = t_s * 1e9 - t0_ns
    r = max(rise_ns, 0.0)
    if t <= 0.0 or t >= flat_ns + 2.0 * r:
        return 0.0
    if r > 0.0 and t < r:
        return t / r
    if t > r + flat_ns:
        return (flat_ns + 2.0 * r - t) / r
    return 1.0


class Device:
    """One configured device; ``run()`` integrates a pulse."""

    def __init__(self, params: ParameterSet = None, mixture: dict = None,
                 electron_rates=None, cavity: Cavity = None,
                 beam_yields: B.BeamYields = None, axis_edges=None):
        p = params if params is not None else engine_parameters()
        self.p = p
        self.mode = p["mode"]
        self.mixture = dict(mixture if mixture is not None
                            else DEFAULT_MIXTURES[self.mode])
        # validate the mixture HERE (the PIC_MCC_3D GasMixture rules): a
        # typo'd gas name must raise, not silently vanish from the fill,
        # and a fraction set that is not a mixture must be refused.
        from .species import REGISTRY
        for g, f in self.mixture.items():
            if g not in REGISTRY or REGISTRY[g].e_form_eV != 0.0:
                raise KeyError(f"mixture gas {g!r} is not a ground-state "
                               "species in the registry")
            if f < 0.0:
                raise ValueError(f"mixture fraction for {g!r} is negative")
        tot = sum(self.mixture.values())
        if not 0.5 < tot < 1.5:
            raise ValueError(
                f"mole fractions sum to {tot:.3f}, which is not a mixture "
                "— give fractions, not percentages or densities")
        self.n_tot = gas_density_m3(p["pressure_atm"], p["T_gas_K"])
        # --- electron-impact tier ------------------------------------------
        if electron_rates is None:
            if p["electron_tier"] == "kinetic":
                from .kinetic_tier import KineticRates
                electron_rates = KineticRates.load_default()
            elif p["electron_tier"] == "maxwell":
                # the Maxwellian-EEDF VERIFICATION tier (0.1.37): Te-indexed
                # only, so the circuit / E/N-direct lane is refused here
                # by capability, with the fix stated (rate-level
                # comparison for the discharge: tools/maxwell_verify.py)
                if (p["drive"] == "discharge"
                        and p["discharge_mode"] == "circuit"):
                    raise ValueError(
                        "the MAXWELL tier is Te-indexed only (no field -> "
                        "EEDF relation) and cannot drive the circuit lane; "
                        "run the circuit on the LOKI tier and compare the "
                        "discharge at rate level with tools/maxwell_verify.py")
                from .loki_tier import rates_for_mode
                from .maxwell_tier import MaxwellShaped
                electron_rates = MaxwellShaped(
                    rates_for_mode(self.mode, self.mixture,
                                   subfloor=bool(p["loki_subfloor"])),
                    self.mixture)
            else:
                from .loki_tier import rates_for_mode
                electron_rates = rates_for_mode(self.mode, self.mixture,
                                                subfloor=bool(p["loki_subfloor"]))
        # beam-mode sub-table bridge (the L3 closure): below the swarm-table
        # floor the MEASURED relaxing-population curves replace the
        # flat+taper stand-in (see kinetic_tier.RelaxationBridged). Only
        # k(Te) lookups change, only below the floor, only where measured;
        # the circuit drive's E/N-direct path delegates through untouched.
        # (the MAXWELL tier is defined at every Te and takes no bridge —
        # its comparison must isolate the shape assumption alone)
        if (p["drive"] == "beam" and p["beam_subtable"] == "relax"
                and p["electron_tier"] != "maxwell"):
            from .kinetic_tier import KineticRates, RelaxationBridged
            relax_src = (electron_rates
                         if isinstance(electron_rates, KineticRates)
                         else KineticRates.load_default())
            electron_rates = RelaxationBridged(electron_rates, relax_src)
        self.tiers = electron_rates
        # --- network + layout ----------------------------------------------
        self.optics = mode_optics(self.mode)
        self.net = build_network(self.mode, self.tiers.available_buckets(),
                                 chen_errata=bool(p["chen_errata"]))
        names = (ARF_MODE_SPECIES if self.mode.startswith("arf")
                 else KRF_MODE_SPECIES)
        self.shape = (p["nx"], p["ny"], p["nz"])
        self.n_cells = int(np.prod(self.shape))
        self.circuit = None
        n_aux = 0
        if p["drive"] == "discharge" and p["discharge_mode"] == "circuit":
            if int(np.prod(self.shape)) != 1:
                raise ValueError("the circuit drive is 0-D in v1: use "
                                 "nx = ny = nz = 1 (the gap is one field)")
            # capability check (0.1.25): the circuit path looks rates and
            # conduction up by E/N — the tables' NATIVE index. The KINETIC
            # tier exports k(<eps>) only (no E/N-indexed curves, no drift
            # velocity in the npz), so this CLI-reachable combination died
            # with an AttributeError on the first RHS call. Refuse it at
            # construction, by capability rather than by name, with the
            # fix stated.
            need = ("k_by_en", "drift_velocity", "en_range_td")
            missing = [m for m in need if not hasattr(self.tiers, m)]
            if missing:
                raise ValueError(
                    "the circuit drive needs E/N-native lookups "
                    f"{need} and the configured electron tier "
                    f"({type(self.tiers).__name__}) lacks {missing} — "
                    "the KINETIC tier exports k(<eps>) only. Run the "
                    "circuit on the LOKI tier (electron_tier='loki').")
            n_aux = N_CIRCUIT_STATES
        self.lay = StateLayout(names, self.n_cells, N_CHANNELS, n_aux)
        self.compiled = CompiledNetwork(
            self.net, self.lay, self.tiers, self.optics,
            petrov17=bool(p["petrov17_f2"]) and p["drive"] == "beam")
        # petrov17 P-normalization needs the fill pressure ([0.1.36])
        self._p17_atm = float(p["pressure_atm"])
        # --- geometry / cavity ---------------------------------------------
        self.extents = np.array([p["Lx_m"], p["Ly_m"], p["Lz_m"]])
        self.d_cell = self.extents / np.array(self.shape, dtype=float)
        self.cell_vol = float(np.prod(self.d_cell))
        self.cell_vols = np.full(self.n_cells, self.cell_vol)
        # --- MESH.txt rung: a refined (nonuniform) beam-axis grid --------
        self.axis_widths = None
        if axis_edges is not None:
            if p["drive"] != "beam":
                raise ValueError("axis_edges refines the BEAM axis — beam "
                                 "drive only")
            if bool(p["transport"]):
                raise ValueError("transport assumes uniform spacing in v1 "
                                 "— not combinable with axis_edges")
            ax = int(p["beam_axis"])
            edges = np.asarray(axis_edges, dtype=float)
            if edges.size != self.shape[ax] + 1 or edges[0] != 0.0 or \
                    abs(edges[-1] - self.extents[ax]) > 1e-12:
                raise ValueError("axis_edges must span [0, L_axis] with "
                                 "n_axis+1 monotone values")
            w = np.diff(edges)
            if (w <= 0).any():
                raise ValueError("axis_edges must be strictly monotone")
            self.axis_widths = w
            trans = self.cell_vol / float(self.d_cell[ax])   # transverse area
            grid = np.zeros(self.shape)
            shape1 = [1, 1, 1]
            shape1[ax] = self.shape[ax]
            grid += w.reshape(shape1) * trans
            self.cell_vols = grid.reshape(-1)
        area = self.extents[0] * self.extents[1]
        if cavity is None:
            # resonator design params (0 / negative sentinels keep the
            # mode defaults — the r_oc/l_op declarations in params.py)
            r_oc = float(p["r_oc"])
            l_op = float(p["l_op"])
            if self.mode.startswith("arf"):
                cavity = Cavity(l_c_m=1.05, l_g_m=self.extents[2],
                                K_axial=1e-5,
                                R_OC=(r_oc if r_oc > 0.0 else 0.08),
                                k_c_loss=80.0,
                                L_op_fixed=(l_op if l_op >= 0.0 else 0.0),
                                area_m2=area, wavelength_nm=193.0,
                                sigma_se_m2=self.optics["sigma_se_m2"],
                                absorbers=dict(self.optics["absorbers"]))
            else:
                cavity = Cavity(l_c_m=self.extents[2], l_g_m=self.extents[2],
                                K_axial=1e-5,
                                R_OC=(r_oc if r_oc > 0.0 else 0.1),
                                k_c_loss=0.0,
                                L_op_fixed=(l_op if l_op >= 0.0 else 0.088),
                                area_m2=area,
                                wavelength_nm=248.0,
                                sigma_se_m2=self.optics["sigma_se_m2"],
                                absorbers=dict(self.optics["absorbers"]))
        self.cav = cavity
        v_g = self.cell_vols.sum()
        if abs(v_g - self.cav.gain_vol_m3) > 1e-9 * v_g:
            raise ValueError(
                f"grid volume {v_g:.6e} != cavity gain volume "
                f"{self.cav.gain_vol_m3:.6e} — the photon-energy bookkeeping "
                "requires them equal (build the cavity from the grid)")
        self.hv_eV = photon_energy_J(self.cav.wavelength_nm) / EV
        self.compiled.set_photon_energy(self.hv_eV)
        if n_aux:
            # gap field along x (Lx = Chen's 12 mm gap); area = Ly x Lz
            self.circuit = CircuitSpec(
                c0_F=p["c0_nF"] * 1e-9, cp_F=p["cp_nF"] * 1e-9,
                lm_H=p["lm_uH"] * 1e-6, ld_H=p["ld_nH"] * 1e-9,
                v0_V=p["v0_kV"] * 1e3, r_ohm=p["r_loop_ohm"],
                gap_m=float(self.extents[0]),
                area_m2=float(self.extents[1] * self.extents[2]),
                rf_ohm=p["rf_ohm"], lf_H=p["lf_nH"] * 1e-9,
                cf_F=p["cf_pF"] * 1e-12)
        # --- drive ----------------------------------------------------------
        self.spec = B.BeamSpec(voltage_kV=p["beam_kV"],
                               j_A_cm2=p["beam_j_A_cm2"],
                               t_start_ns=p["beam_t0_ns"],
                               t_flat_ns=p["beam_flat_ns"],
                               t_rise_ns=p["beam_rise_ns"],
                               axis=int(p["beam_axis"]))
        # optional second identical pulse (rep-rate/inter-pulse studies):
        # same spec time-shifted, so deposition profile and yields are
        # shared; the envelope sum lands in the same drive_in ledger row
        self._beam2_shift_s = 0.0
        if float(p["beam2_t0_ns"]) > 0.0:
            t1_end = (p["beam_t0_ns"] + p["beam_flat_ns"]
                      + 2.0 * p["beam_rise_ns"])
            if p["beam2_t0_ns"] < t1_end:
                raise ValueError(
                    f"beam2_t0_ns={p['beam2_t0_ns']:g} ns overlaps pulse 1 "
                    f"(ends {t1_end:g} ns) — the second pulse must clear "
                    "the first")
            self._beam2_shift_s = (p["beam2_t0_ns"] - p["beam_t0_ns"]) * 1e-9
        self.beam_yields = None
        self._dep_frac = None
        if p["drive"] == "beam":
            dens = {g: f * self.n_tot for g, f in self.mixture.items()}
            if beam_yields is not None:
                self.beam_yields = beam_yields
            else:
                # gases lacking an LXCat set (the 1e-4 Xe trace) still stop
                # the beam but are excluded from the yield split — a <=0.1 %
                # renormalization, stated here
                cs_dens = {g: n for g, n in dens.items() if _has_lxcat(g)}
                try:
                    self.beam_yields = B.spencer_fano_yields(cs_dens)
                except (FileNotFoundError, ValueError) as e:
                    self.beam_yields = B.w_table_yields(self.mixture)
                    self.beam_yields.cite += f" [SF route unavailable: {e}]"
            ax = self.spec.axis
            n_ax = self.shape[ax]
            cells_along = [dens] * n_ax
            dx_axis = (self.axis_widths if self.axis_widths is not None
                       else float(self.d_cell[ax]))
            frac_axis = B.deposition_profile(self.spec, cells_along,
                                             dx_axis)
            # extrude the axial profile over the grid -> per-cell fraction
            grid = np.zeros(self.shape)
            sl = [None, None, None]
            sl[ax] = slice(None)
            shape1 = [1, 1, 1]
            shape1[ax] = n_ax
            grid += frac_axis.reshape(shape1)
            # each transverse column carries frac/n_transverse of the energy
            n_trans = self.n_cells // n_ax
            self._dep_frac = grid.reshape(-1) / n_trans
            self._beam_area = float(np.prod(
                [self.extents[a] for a in range(3) if a != ax]))
        # --- microwave reheat window (a beam-drive overlay; params mw_*) --
        # Physics: during the window a prescribed INJECTED power feeds the
        # overmoded cavity (uniform-field idealization). The field is not
        # a free knob — each step it settles at the lowest E/N where the
        # swarm power balance absorbs the feed,
        #     P_abs = e E v_d(E/N) x nu_m^2/(nu_m^2+omega^2) x N_e,tot ,
        # capped at mw_en_max_td (the empty-cavity ceiling the wall Q
        # supports; a plasma too thin to absorb leaves the excess
        # reflected — only ABSORBED power enters du_e and drive_in).
        # This is the circuit drive's swarm-conduction Joule form with the
        # generator replaced by a power source: dense plasma -> field
        # sags, per-electron heating vanishes; dilute plasma -> field
        # rides the ceiling and absorption dies with n_e. Collision rates
        # stay on the drive's <eps> = (3/2)Te keying: heating drives <eps>
        # toward the swarm point of the settled E/N, where the in-table
        # swarm rates are the field-driven-EEDF physics.
        self._mw_on = float(p["mw_power_W"]) > 0.0
        self._mw_omega = 2.0 * math.pi * float(p["mw_freq_GHz"]) * 1e9
        self._mw_cell_w = None       # stated deposition weights ([0.1.34])
        if self._mw_on and p["drive"] != "beam":
            raise ValueError(
                "the microwave reheat window is a beam-drive overlay "
                "(drive='beam'); the discharge drives already own the "
                "field pathway")
        self._mw_swarm = None
        if self._mw_on:
            # swarm TRANSPORT for the window (v_d(E/N), <eps>(E/N), the
            # E/N range) is an E/N-native quantity separate from the
            # electron-impact RATE tier: the LOKI mixture tables carry it
            # whatever tier keys the rates (the KINETIC tier exports
            # rates only) — stated, and what lets the window run on
            # either tier for the cross-tier check
            if hasattr(self.tiers, "drift_velocity"):
                self._mw_swarm = self.tiers
            else:
                from .loki_tier import rates_for_mode
                self._mw_swarm = rates_for_mode(self.mode, self.mixture)
        if self._mw_on:
            # g(E/N) = (E/N) v_d fac [Td m/s] on a log grid up to the
            # ceiling; the rhs picks the LOWEST crossing (the field
            # builds from below — the stable branch through any NDC
            # region of the Ar drift curve)
            en_max = float(p["mw_en_max_td"])
            grid = np.geomspace(1e-3, en_max, 512)
            g = np.empty_like(grid)
            for i, en in enumerate(grid):
                v_d = self._mw_swarm.drift_velocity(float(en))
                e_f = en * 1e-21 * self.n_tot
                if v_d <= 0.0:
                    g[i] = 0.0
                    continue
                nu = E_CHARGE_SIGNED * e_f / (M_E * v_d)
                g[i] = en * v_d * nu * nu / (nu * nu + self._mw_omega ** 2)
            self._mw_grid_en = grid
            self._mw_grid_g = g
            self._mw_g_ceiling = float(g[-1])
        # --- in-loop ASE precomputation (S6): centre-point chords over a
        # coarse Gauss-Legendre x phi quadrature, geometry-static; per RHS
        # call only the closed-form uniform-medium intensity is evaluated
        # on the VOLUME-AVERAGED g0/alpha (stated approximations: averaged
        # field, centre-point flux; the axial cavity solid angle ~K_axial
        # overlaps the quadrature — a 1e-5-class double count, stated)
        self._ase_on = bool(p["ase_inloop"])
        # lazily built on the first `jac` call (0.1.39)
        self._jac_struct = None
        self._jac_groups = None
        self._rhs_cache = None
        self._jac_factor = _UNSET
        self._jac_factor_cap = 1e12    # see `jac` (0.1.39)
        self._jac_floor = False        # opt-in; see `jac` (0.1.39)
        from .ase import sphere_quadrature, sum_rule_bandwidth_hz, \
            planck_source
        dirs, wts = sphere_quadrature(6, 8)
        c = 0.5 * self.extents
        with np.errstate(divide="ignore", invalid="ignore"):
            t_hi = np.where(dirs > 0, (self.extents - c) / dirs, np.inf)
            t_lo = np.where(dirs < 0, -c / dirs, np.inf)
        self._ase_L = np.min(np.where(np.isnan(np.minimum(t_hi, t_lo)),
                                      np.inf, np.minimum(t_hi, t_lo)),
                             axis=1)
        self._ase_w = wts
        self._ase_sbd = (planck_source(self.cav.wavelength_nm)
                         * sum_rule_bandwidth_hz(self.cav.wavelength_nm,
                                                 self.cav.sigma_se_m2,
                                                 self.optics["tau_upper_s"]))
        self._abs_idx = [(self.lay.index[sp], sig)
                         for sp, sig in self.cav.absorbers.items()
                         if sp in self.lay.index]
        # --- transport stencil (optional; conservative 6-neighbour) -------
        self._transport = bool(p["transport"])
        if self._transport:
            self._diff0 = (p["diff_m2_s"] * (p["T_gas_K"] / 300.0) ** 1.5
                           / p["pressure_atm"])
        # per-species production per pair, mapped to layout rows once — and
        # MASS-CONSERVING: each produced species consumes its ground state
        # (Ar+ eats an Ar, F2(a3) eats an F2; 0.1.25: k product units may
        # consume ONE ground of k-fold composition — 2 F eat one F2, the
        # prompt-dissociation booking). Found mechanically: the layout
        # species with the same elemental atoms (or the atoms x k), charge
        # 0 and E_form 0.
        self._beam_rows = []
        if self.beam_yields is not None:
            from .species import atoms_of, charge_of, e_form_of
            grounds = {}
            for s in self.lay.species:
                if charge_of(s) == 0 and e_form_of(s) == 0.0:
                    grounds[tuple(sorted(atoms_of(s).items()))] = s
            for s, n_pp in self.beam_yields.per_pair.items():
                i = self.lay.index.get(s)
                if i is None or n_pp <= 0.0:
                    continue
                self._beam_rows.append((i, n_pp))
                atoms = atoms_of(s)
                g = grounds.get(tuple(sorted(atoms.items())))
                mult = 1
                if g is None:
                    for k in (2, 3):
                        scaled = tuple(sorted((el, k * v)
                                              for el, v in atoms.items()))
                        g = grounds.get(scaled)
                        if g is not None:
                            mult = k
                            break
                if g is None:
                    raise KeyError(f"beam product {s!r} has no ground-state "
                                   "species in the layout to consume")
                self._beam_rows.append((self.lay.index[g], -n_pp / mult))

    def _mw_field_td(self, ne_weighted: float, p_inj_W: float) -> float:
        """The settled window E/N [Td]: lowest grid crossing of the swarm
        power balance, else the ceiling. ne_weighted = sum(n_e V) [m^-3
        m^3 = electrons]."""
        if p_inj_W <= 0.0:
            return 0.0
        c = E_CHARGE_SIGNED * 1e-21 * self.n_tot * max(ne_weighted, 0.0)
        if c <= 0.0:
            return float(self._mw_grid_en[-1])
        g_need = p_inj_W / c
        if g_need >= self._mw_g_ceiling:
            return float(self._mw_grid_en[-1])
        # first crossing (argmax on the mask: g need not be monotone
        # through the Ar NDC dip — the lowest crossing is the branch the
        # field reaches building up from zero)
        j = int(np.argmax(self._mw_grid_g >= g_need))
        if j <= 0:
            return float(self._mw_grid_en[0])
        g0, g1 = self._mw_grid_g[j - 1], self._mw_grid_g[j]
        e0, e1 = self._mw_grid_en[j - 1], self._mw_grid_en[j]
        w = (g_need - g0) / (g1 - g0) if g1 > g0 else 0.0
        return float(e0 + w * (e1 - e0))

    def set_mw_cell_weights(self, w) -> None:
        """Stated per-cell microwave DEPOSITION weights (a skin/shell
        profile, or a measured / full-wave deposition map). The window's
        power balance settles against the weighted conductivity and
        deposits per cell as w_i n_e_i, so the PRESCRIBED injected power
        is honored to the settle grid's interpolation accuracy (~1e-6
        relative) for any non-negative, non-zero profile — the settled
        field rescales, and the ledger books the ACTUAL deposit exactly
        (gate: test_mw_cell_weights_identity_masking_and_conservation).
        Default (None) = conductivity weighting at the common settled
        field — the uniform-field limit, bit-identical to the pre-hook
        behavior. [0.1.34]"""
        if w is None:
            self._mw_cell_w = None
            return
        w = np.asarray(w, dtype=float).reshape(-1)
        if w.shape != (self.n_cells,) or (w < 0.0).any() or w.sum() <= 0.0:
            raise ValueError("mw cell weights must be one non-negative "
                             "value per cell with a positive sum")
        self._mw_cell_w = w

    # ------------------------------------------------------------------ state
    def initial_state(self) -> np.ndarray:
        p = self.p
        y = np.zeros(self.lay.size)
        cells = self.lay.cells_view(y)
        dens = cells[: self.lay.n_species]
        dens[:] = p["seed_density_m3"]
        for g, frac in self.mixture.items():
            i = self.lay.index.get(g)
            if i is not None:
                dens[i] = frac * self.n_tot
        if p["drive"] == "discharge":
            buffer_ion = "Ne+" if self.mixture.get("Ne", 0) >= \
                self.mixture.get("Ar", 0) else "Ar+"
            i = self.lay.index[buffer_ion]
            dens[i] = np.maximum(dens[i], p["preion_m3"])
        if self.circuit is not None:
            y[self.lay.i_aux0:self.lay.i_aux0 + N_CIRCUIT_STATES] = \
                self.circuit.initial()
        ne0 = self.lay.electron_density(dens)
        cells[self.lay.i_ue] = 1.5 * ne0 * p["te_seed_eV"] * EV
        cells[self.lay.i_tg] = p["T_gas_K"]
        y[self.lay.i_ph] = 1.0
        return y

    def te_of(self, y) -> np.ndarray:
        cells = self.lay.cells_view(y)
        dens = np.maximum(cells[: self.lay.n_species], 0.0)
        # 1 m^-3 floor: physically nothing, numerically decisive — the
        # undriven afterglow past ~500 ns reaches n_e -> 0 and the bare
        # 0/0 NaN'd the numerical Jacobian (found by the 0.1.21 keeper
        # runs; the clip cannot rescue a NaN)
        ne = np.maximum(self.lay.electron_density(dens), 1.0)
        return np.clip((2.0 / 3.0) * cells[self.lay.i_ue] / (ne * EV),
                       0.026, 50.0)

    # -------------------------------------------------------------------- rhs
    def rhs(self, t, y):
        lay = self.lay
        p = self.p
        dy = np.zeros_like(y)
        cells = lay.cells_view(y)
        d_cells = lay.cells_view(dy)
        dens = np.maximum(cells[: lay.n_species], 0.0)
        # same 1 m^-3 floor as te_of (deep-afterglow 0/0 guard)
        ne = np.maximum(lay.electron_density(dens), 1.0)
        te = np.clip((2.0 / 3.0) * cells[lay.i_ue] / (ne * EV), 0.026, 50.0)
        tg = np.clip(cells[lay.i_tg], 200.0, 3000.0)
        n_ph = float(y[lay.i_ph])
        en_td = None
        E_gap = 0.0
        if self.circuit is not None:
            E_gap = float(y[lay.i_aux0 + IE])
            en_td = abs(E_gap) / self.n_tot / 1e-21
        ase_rate = self._ase_rate(dens) if self._ase_on else 0.0
        # beam deposition is computed BEFORE the chemistry call so the
        # petrov17 beam-EEDF rates ([0.1.36], flag-gated) can see the
        # instantaneous per-cell power; the drive block below reuses the
        # same arrays unchanged
        p_dep = None
        if p["drive"] == "beam":
            p_areal = self.spec.power_areal_W_m2(t)
            if self._beam2_shift_s > 0.0:
                p_areal += self.spec.power_areal_W_m2(t - self._beam2_shift_s)
            if p_areal > 0.0:
                e_rate = p_areal * self._beam_area          # W into the grid
                p_dep = e_rate * self._dep_frac / self.cell_vols  # W/m^3
        p17_norm = None
        if p_dep is not None and self.compiled._petrov17:
            # Petrov 2017: P = 2 P_beam[MW/cm^3]/p[atm] + 1e-3
            p17_norm = 2.0 * (p_dep / 1.0e12) / self._p17_atm + 1.0e-3
        ch = self.compiled.rhs(dens, te, tg, ne, n_ph, self.n_tot,
                               en_td=en_td, ase_rate=ase_rate,
                               p17_norm=p17_norm)
        ddens = ch["ddens"]
        du_e = ch["du_e"]
        dgas = ch["dgas"]
        drive_rate = 0.0
        # --- drive ----------------------------------------------------------
        if p["drive"] == "beam":
            if p_dep is not None:
                pairs = p_dep / (self.beam_yields.W_eV * EV)
                for i_s, n_pp in self._beam_rows:
                    ddens[i_s] += n_pp * pairs
                du_e += pairs * self.beam_yields.eps_sec_eV * EV
                dgas += pairs * self.beam_yields.subfloor_eV * EV
                drive_rate = float(p_dep @ self.cell_vols)
        elif self.circuit is not None:
            # E/N-driven Joule heating: j = e n_e v_d(E/N) sign(E); the
            # electron gas receives j.E >= 0 and drive_in integrates
            # exactly that (the circuit's stores are upstream of the cell
            # budget — circuit.py docstring)
            cs = self.circuit
            v_d = self.tiers.drift_velocity(en_td)
            j_cond = (E_CHARGE_SIGNED * float(ne[0]) * v_d
                      * math.copysign(1.0, E_gap))
            p_joule = j_cond * E_gap                     # W/m^3, >= 0
            du_e += p_joule
            drive_rate = p_joule * float(self.cell_vols.sum())
            d_aux = cs.derivs(y[lay.i_aux0:lay.i_aux0 + N_CIRCUIT_STATES],
                              j_cond * cs.area_m2)
            dy[lay.i_aux0:lay.i_aux0 + N_CIRCUIT_STATES] = d_aux
        else:
            x = (t * 1e9 - p["dis_t0_ns"]) / max(p["dis_fwhm_ns"], 1e-3)
            p_dis = p["dis_peak_W_m3"] * math.exp(-2.772588722 * x * x)
            du_e += p_dis
            drive_rate = p_dis * float(self.cell_vols.sum())
        # --- microwave reheat window (beam overlay; init block above) -----
        if self._mw_on:
            env = _mw_envelope(t, p["mw_t0_ns"], p["mw_rise_ns"],
                               p["mw_flat_ns"])
            if env > 0.0:
                p_inj = float(p["mw_power_W"]) * env
                ne_eff = (ne if self._mw_cell_w is None
                          else ne * self._mw_cell_w)
                en_mw = self._mw_field_td(float(ne_eff @ self.cell_vols),
                                          p_inj)
                v_d = self._mw_swarm.drift_velocity(en_mw)
                if v_d > 0.0:
                    e_mw = en_mw * 1e-21 * self.n_tot        # V/m, RMS
                    nu_m = E_CHARGE_SIGNED * e_mw / (M_E * v_d)
                    fac = nu_m * nu_m / (nu_m * nu_m + self._mw_omega ** 2)
                    p_mw = E_CHARGE_SIGNED * ne_eff * v_d * e_mw * fac  # W/m^3
                    du_e += p_mw
                    drive_rate += float(p_mw @ self.cell_vols)
        # --- photon equation (gain-volume averages) ------------------------
        w = self.cell_vols
        v_g = float(w.sum())
        r_ste_avg = float(ch["r_ste"] @ w) / v_g
        r_x_avg = float(ch["r_x"] @ w) / v_g
        r_se_avg = float(ch["r_se"] @ w) / v_g
        i_f2 = lay.index.get("F2")
        c_f2 = (float(dens[i_f2] @ w) / v_g / self.n_tot
                if i_f2 is not None else 0.0)
        dy[lay.i_ph] = photon_rhs(self.cav, r_ste_avg, r_x_avg, r_se_avg,
                                  n_ph, c_f2)
        # --- ledger channels ------------------------------------------------
        hv_J = self.cav.photon_J
        v_cav = self.cav.cavity_vol_m3
        out_rate, int_rate = loss_rates(self.cav, n_ph, c_f2)
        spont_total = float(ch["r_se"] @ w) * hv_J           # W emitted
        spont_captured = self.cav.K_axial * r_se_avg * v_cav * hv_J
        off = lay.i_ledger0
        dy[off + 0] = drive_rate
        dy[off + 1] = out_rate * hv_J * v_cav
        dy[off + 2] = int_rate * hv_J * v_cav
        dy[off + 3] = (float(ch["rad_escape"] @ w)
                       + spont_total - spont_captured)
        if self._transport and self.n_cells > 1:
            ddens = ddens + self._transport_term(dens)
        # --- assemble cell derivatives -------------------------------------
        d_cells[: lay.n_species] = ddens
        d_cells[lay.i_ue] = du_e
        d_cells[lay.i_tg] = dgas / (1.5 * self.n_tot * K_B)
        # --- non-finite guard (0.1.25): an inf/NaN derivative is never
        # legitimate — without this it surfaced as an inscrutable
        # LinAlgError deep in the BDF lu_factor (found running the Chen
        # AS-PRINTED rate table under the paper's true base-case circuit).
        # Raise with the state location instead, so a blow-up names its
        # physics. One isfinite pass per call; cost negligible.
        if not np.isfinite(dy).all():
            bad = np.flatnonzero(~np.isfinite(dy))
            names = []
            for b in bad[:4]:
                if b < lay.per_cell * lay.n_cells:
                    r = int(b % lay.per_cell)
                    names.append(lay.species[r] if r < lay.n_species else
                                 ("U_e" if r == lay.i_ue else "T_g"))
                elif b == lay.i_ph:
                    names.append("n_ph")
                elif b >= lay.i_aux0:
                    names.append(f"circuit[{int(b - lay.i_aux0)}]")
                else:
                    names.append(f"ledger[{int(b - lay.i_ledger0)}]")
            raise FloatingPointError(
                f"non-finite RHS at t = {t * 1e9:.3f} ns in {names} "
                f"({bad.size} entries) — the state has blown up; check the "
                "rate table / drive for a physically impossible input "
                "(e.g. chen_errata=False under the true Chen circuit)")
        # the solver evaluates f at (t, y) immediately before asking for a
        # Jacobian there; caching it saves `jac` a full extra RHS call (the
        # guard is exact, so a stale entry cannot be used)
        self._rhs_cache = (t, y.copy(), dy)
        return dy

    def _ase_rate(self, dens) -> float:
        """Per-upper-molecule ASE depopulation rate [1/s] from the 4-pi
        flux of the volume-averaged gain field (closed-form chords)."""
        w = self.cell_vols
        v_g = float(w.sum())
        i_up = self.lay.index[self.optics["upper"]]
        g0 = self.cav.sigma_se_m2 * float(dens[i_up] @ w) / v_g
        if g0 <= 0.0:
            return 0.0
        alpha = 0.0
        for i_sp, sig in self._abs_idx:
            alpha += sig * float(dens[i_sp] @ w) / v_g
        gn = g0 - alpha
        L = self._ase_L
        x = np.clip(gn * L, -60.0, 60.0)
        small = np.abs(x) < 1e-8
        e1 = np.where(small, L * (1.0 + 0.5 * x),
                      np.expm1(x) / np.where(small, 1.0, gn))
        phi = float(np.sum(self._ase_w * g0 * self._ase_sbd * e1))
        return self.cav.sigma_se_m2 * phi / self.cav.photon_J

    def _transport_term(self, dens):
        """Conservative 6-neighbour diffusion divergence [m^-3 s^-1]:
        D (n_j - n_i) A_f / (d_f V_i) summed over faces — exact pairwise
        antisymmetry makes global conservation machine-exact."""
        nx, ny, nz = self.shape
        f = dens.reshape(dens.shape[0], nx, ny, nz)
        out = np.zeros_like(f)
        D = self._diff0
        for ax, h in ((1, self.d_cell[0]), (2, self.d_cell[1]),
                      (3, self.d_cell[2])):
            if f.shape[ax] < 2:
                continue
            d = np.diff(f, axis=ax) * (D / (h * h))   # flux/h per face
            pad = [slice(None)] * 4
            padl = list(pad); padl[ax] = slice(0, -1)
            padr = list(pad); padr[ax] = slice(1, None)
            out[tuple(padl)] += d
            out[tuple(padr)] -= d
        return out.reshape(dens.shape)

    def atol_vector(self):
        """Per-state absolute tolerance. One definition, used both by the
        integrator's error test and as the floor on the Jacobian's
        differencing step (0.1.39) — the two must agree, or the Jacobian
        would resolve changes the error test calls noise."""
        atol = np.full(self.lay.size, 1e3)          # densities [m^-3]
        cells_atol = self.lay.cells_view(atol)
        cells_atol[self.lay.i_ue] = 1e-10           # U_e [J/m^3]
        cells_atol[self.lay.i_tg] = 1e-6            # T_g [K]
        atol[self.lay.i_ph] = 1e-3
        atol[self.lay.i_ledger0:] = 1e-12
        if self.lay.n_aux:
            # q0, qp [C], i1, i2 [A], E [V/m], if [A], qf [C]
            atol[self.lay.i_aux0:self.lay.i_aux0 + N_CIRCUIT_STATES] = \
                [1e-12, 1e-12, 1e-6, 1e-6, 1e-2, 1e-6, 1e-15]
        return atol

    # ------------------------------------------------------------- Jacobian
    def jac(self, t, y):
        """Jacobian for the BDF Newton matrix: scipy's adaptive ``num_jac``
        with the structurally dead columns pinned and the step factor capped.
        (The first 0.1.39 draft used a FIXED step; it was measured 7.4x
        slower and replaced — see the comment block below. Merged into the
        Agent tree in 0.1.40.)

        0.1.39. Wraps scipy's ``num_jac``, which carries an adaptive
        per-column step ``factor`` and multiplies it by 10 every time a
        column's finite difference comes back "too small". A column that is
        structurally dead — every ledger column here, since no rate anywhere
        reads an accumulated ledger value — is too small on EVERY call, so
        its factor grows without bound and eventually overflows against a
        density scale of ~1e25 (``scipy/integrate/_ivp/common.py`` in
        ``_sparse_num_jac``: ``new_factor * y_scale[ind]``). That is the
        0.1.35 stiff blocker: it killed the 0.8 % F2 lever point and the
        multi-cell 3-D runs, and it is a solver-bookkeeping fault, not a
        physics one. Handing BDF a Jacobian means ``num_jac`` never runs.

        The derivative itself is the same first-order difference BDF would
        have taken, evaluated through ``scipy.optimize._numdiff`` with the
        textbook ``sqrt(eps)`` relative step and the engine's own sparsity
        groups, so the cost is one RHS call per colour group exactly as
        before. Nothing about the equations changes; only the machinery that
        differentiates them, and a Newton matrix does not enter the error
        test (the step controller still governs accuracy).
        """
        from scipy.integrate._ivp.common import num_jac
        if self._jac_factor is _UNSET:
            from scipy.optimize._numdiff import group_columns
            from scipy.sparse import csc_matrix
            struct = csc_matrix(self.jac_sparsity())
            # Sparsity goes to num_jac only for a real grid, exactly as
            # before 0.1.39. Handing it a 0-D system puts BDF on its SPARSE
            # linear algebra for a ~42x42 matrix whose photon row is dense,
            # and sparse LU on that is far slower than the dense path it
            # used to take: measured, test_engine.py went from part of a
            # 12-minute suite to over 10 minutes on its own.
            self._jac_struct = ((struct, group_columns(struct))
                                if self.n_cells > 1 else None)
            self._jac_atol = self.atol_vector()
            self._jac_factor = None
            # the columns nothing depends on: exactly the pathology below
            self._jac_dead = np.asarray(
                np.diff(struct.indptr) == 0).ravel()
        c = self._rhs_cache
        f0 = (c[2] if (c is not None and c[0] == t
                       and np.array_equal(c[1], y)) else self.rhs(t, y))
        J, factor = num_jac(self._rhs_vectorized, t, y, f0,
                            self._jac_atol, self._jac_factor,
                            self._jac_struct)
        # THE CAP is the whole fix. num_jac's per-column `factor` is the
        # relative step; it multiplies by ten whenever a column's difference
        # comes back "too small" and there is NO UPPER BOUND, so a
        # structurally dead column climbs until `factor * y_scale` overflows
        # against a ~1e25 density. The cap is therefore set FAR out, at 1e12:
        # `factor * y_scale` then tops out near 1e37, ~270 decades below the
        # overflow, while every legitimate growth a live column asks for is
        # still granted. A tight cap is NOT safe here — capping at 1 looked
        # fine on a 40 ns multi-cell run (parity, 6 Jacobians) and cost 12x on
        # the same case run to 100 ns (test_3d_gradient: 4259 s against 342 s),
        # because in the AFTERGLOW collapsing densities legitimately need
        # steps far above their own magnitude and the cap was refusing them.
        # Pinning the dead columns is what stops the runaway; the cap only
        # backstops it.
        #
        # Both halves of that sentence are measured, on the 24-cell gradient
        # case against stock scipy (7.8 s, 176 steps, 6 Jacobians):
        #   fixed sqrt(eps) step, no adaptation : 70.9 s, 296 steps, 50 jac
        #   cap 1e-2 (a 1 % step ceiling)       : 13.1 s, 181 steps,  9 jac
        #   cap 1e0  (40-ns window; see above)  :  8.5 s, 176 steps,  6 jac
        # (the shipped cap is 1e12: 1e0 was then measured 12x slow at 100 ns)
        # Killing the adaptation to kill the runaway costs 7.4x. Bounding it
        # costs nothing. The first version of this fix made that mistake and
        # took the suite from 12 minutes to 5 hours before it was measured.
        # FLOOR as well as cap. num_jac's step is h = factor * max(|y|, atol);
        # forcing h >= atol keeps a near-zero column's difference out of the
        # rounding noise, which is what the near-threshold cavity cases need,
        # while the cap keeps the runaway impossible. Both bounds are inert
        # for a well-scaled column, where the adaptation runs untouched.
        # Pin the DEAD columns, leave every live one alone. This is the whole
        # fix and it is deliberately surgical: a column with no structural
        # nonzeros is the only one that can never produce a difference, so it
        # is the only one whose factor climbs on every single call. Live
        # columns self-correct — their difference stops being "too small" as
        # soon as the step is big enough — and their adaptation is what makes
        # the algorithm work at all, so it is left exactly as scipy wrote it.
        # The cap is belt-and-braces on top.
        factor = np.clip(factor, 1e-14, self._jac_factor_cap)
        factor[self._jac_dead] = np.sqrt(np.finfo(float).eps)
        # A state at or below its OWN tolerance is indistinguishable from
        # zero, and a sqrt(eps) step off it lands in the rounding noise, so
        # probe it with a full tolerance instead (y_scale is atol there, so
        # factor 1 IS a step of one atol). The bound is deliberately narrow:
        # applying it to every state instead — including species merely
        # small, at 1e4..1e8 m^-3, where it means a 0.1..10 % perturbation —
        # was measured to cost 12x on the 24-cell case while helping the
        # 0-D one, which is what makes "below its own atol" the right line.
        # TOLERANCE FLOOR — OFF by default, opt in with `dev._jac_floor = True`.
        # Forcing every step to h >= atol is the one thing measured to rescue
        # a near-threshold 0-D cavity run: the 0.8 % F2 lever point goes from
        # not finishing in ten minutes to 8 s, same answer to four figures.
        # It is expensive everywhere else and must not be a default — 12x on
        # a 24-cell grid (121 s against 9.8 s), and 6x across the whole test
        # suite when it was briefly defaulted on for 0-D. The reason is that
        # a tolerance-scale probe is a large RELATIVE perturbation for any
        # species that is merely small, and it degrades the matrix; it pays
        # only where a near-threshold cavity field makes the alternative a
        # step lost in rounding noise. Narrower forms do not work: flooring
        # only states already below their atol, and flooring only the photon
        # column, both leave that case unfinished.
        if self._jac_floor:
            np.maximum(factor,
                       self._jac_atol / np.maximum(np.abs(y), self._jac_atol),
                       out=factor)
        self._jac_factor = factor
        return J

    def _rhs_vectorized(self, t, z):
        """`num_jac` probes a whole colour group in one call: z is (n, k)."""
        if z.ndim == 1:
            return self.rhs(t, z)
        out = np.empty_like(z)
        for i in range(z.shape[1]):
            out[:, i] = self.rhs(t, z[:, i])
        return out

    # ------------------------------------------------------------- sparsity
    def jac_sparsity(self):
        """Block-diagonal cell pattern + photon/aux coupling — PLUS the two
        GLOBAL volume-average couplings the 0.1.24 pattern missed (found
        2026-08-25): the in-loop ASE rate (volume-averaged upper/absorber
        field -> every cell's upper/fragment/T_g rows) and the MW window
        field (volume-summed n_e -> every cell's U_e row). A missing
        sparsity entry starves the BDF Newton matrix of a real derivative —
        convergence quality, not answer correctness (the error test still
        governs accuracy) — and both couplings are cheap supersets here."""
        from scipy import sparse
        n = self.lay.size
        pc = self.lay.per_cell
        C = self.lay.n_cells
        m = sparse.lil_matrix((n, n), dtype=np.int8)
        block = np.ones((pc, pc), dtype=np.int8)
        for c in range(C):
            a = c * pc
            m[a:a + pc, a:a + pc] = block
        i_ph = self.lay.i_ph
        m[i_ph, : i_ph + 1] = 1
        m[: i_ph + 1, i_ph] = 1
        for a in range(self.lay.i_aux0, self.lay.size):
            m[a, :] = 1
            m[:, a] = 1
        if C > 1 and self._ase_on:
            i_up = self.lay.index[self.optics["upper"]]
            cols_c = sorted({i_up} | {i for i, _s in self._abs_idx})
            rows_c = sorted({i_up, self.lay.i_tg}
                            | set(self.compiled._frag_idx))
            cols = np.array([c * pc + j for c in range(C) for j in cols_c])
            rows = np.array([c * pc + j for c in range(C) for j in rows_c])
            m[np.ix_(rows, cols)] = 1
        if C > 1 and self._mw_on:
            chg = sorted(set(self.lay._pos_idx) | set(self.lay._neg_idx))
            cols = np.array([c * pc + j for c in range(C) for j in chg])
            rows = np.array([c * pc + self.lay.i_ue for c in range(C)])
            m[np.ix_(rows, cols)] = 1
        if self._transport and C > 1:
            nx, ny, nz = self.shape
            S = self.lay.n_species
            idx = np.arange(C).reshape(nx, ny, nz)
            for ax in range(3):
                if idx.shape[ax] < 2:
                    continue
                a = np.moveaxis(idx, ax, 0)
                for i, j in ((a[:-1].ravel(), a[1:].ravel()),):
                    for ci, cj in zip(i, j):
                        m[ci * pc:ci * pc + S, cj * pc:cj * pc + S] = 1
                        m[cj * pc:cj * pc + S, ci * pc:ci * pc + S] = 1
        return m.tocsr()

    # ------------------------------------------------------------------- run
    def run(self, n_save: int = 241, dense_cells: bool = None) -> dict:
        """Integrate to t_end_ns and collect the run histories.

        ``dense_cells`` keeps the full state at every save time as
        ``res["y_t"]`` (shape ``(lay.size, n_save)``), so per-cell fields
        can be reconstructed after the run (``tools/visualize_3d.py``) —
        the volume-averaged histories in ``res`` summarize honestly, but
        the 3-D structure lives only in ``y_t``. ``None`` (the default)
        resolves to ``n_cells > 1``: a 3-D run keeps its cells (~8 MB at
        240 cells x 121 saves), a 0-D run stays lean.
        """
        dense = ((self.n_cells > 1) if dense_cells is None
                 else bool(dense_cells))
        from scipy.integrate import solve_ivp
        p = self.p
        t_end = p["t_end_ns"] * 1e-9
        y0 = self.initial_state()
        t_eval = np.linspace(0.0, t_end, int(n_save))
        atol = self.atol_vector()
        # 0.1.39: hand BDF an explicit Jacobian (see `jac`) instead of the
        # sparsity pattern alone. The pattern still shapes the differencing,
        # but scipy's own `num_jac` — whose unbounded per-column step factor
        # is the 0.1.35 stiff blocker — never runs.
        kw = {"jac": self.jac}
        t0 = _time.time()
        sol = solve_ivp(self.rhs, (0.0, t_end), y0, method="BDF",
                        rtol=p["rtol"], atol=atol, t_eval=t_eval,
                        max_step=2e-9, **kw)
        wall_s = _time.time() - t0
        out = self._collect(sol, y0, wall_s, dense)
        return out

    def _collect(self, sol, y0, wall_s, dense_cells: bool = False) -> dict:
        lay = self.lay
        n_t = sol.t.size
        dens_avg = {}
        w = self.cell_vols
        v_g = float(w.sum())
        te_avg = np.zeros(n_t)
        tg_avg = np.zeros(n_t)
        ne_avg = np.zeros(n_t)
        g0 = np.zeros(n_t)
        p_out = np.zeros(n_t)
        closure = []
        for s in lay.species:
            dens_avg[s] = np.zeros(n_t)
        from .laser import output_power_W
        for j in range(n_t):
            yj = sol.y[:, j]
            cells = lay.cells_view(yj)
            dens = np.maximum(cells[: lay.n_species], 0.0)
            for s, i in lay.index.items():
                dens_avg[s][j] = float(dens[i] @ w) / v_g
            ne = lay.electron_density(dens)
            ne_avg[j] = float(ne @ w) / v_g
            te_avg[j] = float(self.te_of(yj) @ w) / v_g
            tg_avg[j] = float(cells[lay.i_tg] @ w) / v_g
            davg = {s: dens_avg[s][j] for s in lay.species}
            g0[j] = gain_m1(self.cav, davg, self.optics["upper"])
            p_out[j] = output_power_W(self.cav, float(yj[lay.i_ph]))
            closure.append(closure_report(
                lay, y0, yj, self.cell_vols, self.n_tot, self.p["T_gas_K"],
                self.cav.cavity_vol_m3, self.cav.photon_J))
        from .ledger import channels_J
        final = channels_J(lay, sol.y[:, -1])
        res = {
            "success": bool(sol.success), "t_ns": sol.t * 1e9,
            "dens_avg": dens_avg, "Te_eV": te_avg, "Tg_K": tg_avg,
            "ne_m3": ne_avg, "n_ph": sol.y[lay.i_ph], "g0_m1": g0,
            "P_out_W": p_out, "closure": closure,
            "final_channels_J": final,
            "eta_intrinsic": (final["laser_out"] / final["drive_in"]
                              if final["drive_in"] > 0 else 0.0),
            "wall_s": wall_s, "y_final": sol.y[:, -1].copy(), "y0": y0,
            "device": self,
        }
        if dense_cells:
            res["y_t"] = sol.y.copy()
        # --- validity gauges (the ESSim.stability analog, 2026-08-04):
        # the stated assumptions, each with the number that says how far
        # this run strayed from it ------------------------------------------
        eps_lo = eps_hi = None
        for xs, _ys in self.tiers.curves.values():
            lo, hi = float(np.exp(xs[0])), float(np.exp(xs[-1]))
            eps_lo = lo if eps_lo is None else min(eps_lo, lo)
            eps_hi = hi if eps_hi is None else max(eps_hi, hi)
        eps_tr = 1.5 * te_avg
        i_f2 = lay.index.get("F2")
        f2_tr = dens_avg["F2"] if i_f2 is not None else None
        res_validity = {
            "ionization_frac_max": float(np.max(ne_avg) / self.n_tot),
            "f2_burn_frac": (float(1.0 - f2_tr.min() / max(f2_tr[0], 1e-300))
                             if f2_tr is not None and f2_tr[0] > 0 else 0.0),
            "te_in_table_frac": float(np.mean((eps_tr >= eps_lo)
                                              & (eps_tr <= eps_hi))),
            "tier_eps_range_eV": (eps_lo, eps_hi),
        }
        # the L3 bridge gauge: fraction of the pulse whose <eps> is covered
        # by table OR measured relaxation curve (beam mode with the bridge)
        cov = getattr(self.tiers, "sub_table_coverage", None)
        if cov is not None:
            cov = cov()
            if cov:
                cov_lo = min(lo for lo, _fl in cov.values())
                res_validity["relax_bridged_buckets"] = sorted(cov)
                res_validity["relax_eps_floor_eV"] = cov_lo
                res_validity["te_in_covered_frac"] = float(
                    np.mean((eps_tr >= cov_lo) & (eps_tr <= eps_hi)))
        if self._mw_on:
            en_cap = float(self.p["mw_en_max_td"])
            v_d0 = self._mw_swarm.drift_velocity(en_cap)
            e0 = en_cap * 1e-21 * self.n_tot
            nu0 = (E_CHARGE_SIGNED * e0 / (M_E * v_d0)) if v_d0 > 0.0 else 0.0
            lo, hi = self._mw_swarm.en_range_td
            res_validity["mw_power_W"] = float(self.p["mw_power_W"])
            res_validity["mw_window_ns"] = (
                float(self.p["mw_t0_ns"]),
                float(self.p["mw_t0_ns"] + 2.0 * self.p["mw_rise_ns"]
                      + self.p["mw_flat_ns"]))
            res_validity["mw_en_max_td"] = en_cap
            res_validity["mw_eps_swarm_at_cap_eV"] = float(
                self._mw_swarm.mean_energy_at_en(en_cap))
            res_validity["mw_nu_over_omega_at_cap"] = (
                nu0 / self._mw_omega if nu0 > 0.0 else 0.0)
            res_validity["mw_cap_in_table"] = bool(lo <= en_cap <= hi)
        if self.circuit is not None:
            a0 = lay.i_aux0
            aux = sol.y[a0:a0 + N_CIRCUIT_STATES]
            en_t = np.abs(aux[IE]) / self.n_tot / 1e-21
            lo, hi = self.tiers.en_range_td
            res["circuit"] = {
                "U_gap_V": aux[IE] * self.circuit.gap_m,
                "I2_A": aux[II2], "I1_A": aux[II1],
                "en_td": en_t,
                "en_above_table_frac": float((en_t > hi).mean()),
                # the gauge was ONE-SIDED through 0.1.24: the breakdown
                # ramp (and everything under the table floor — 5 Td in
                # the Chen set, which ships no low-E/N companion) rides
                # the flat+taper stand-in, and that fraction went
                # unreported while the scorecard said "100% in-table"
                "en_below_table_frac": float((en_t < lo).mean()),
                "en_table_range_td": (lo, hi),
                "circuit_energy_J": np.array(
                    [self.circuit.energy_J(aux[:, j])
                     for j in range(aux.shape[1])]),
            }
        res["validity"] = res_validity
        res["regime"] = self.regime_line(res)
        return res

    def regime_line(self, res) -> str:
        y = self.beam_yields
        route = (f"beam yields via {y.mode}" if y is not None else
                 ("circuit discharge (E/N-direct rates, swarm-drift "
                  "conduction)" if self.circuit is not None
                  else "prescribed discharge"))
        if y is not None:
            br = res["validity"].get("relax_bridged_buckets")
            if br:
                route += (", sub-table rates = MEASURED relax curves ("
                          + ",".join(br) + "; <eps> in-coverage "
                          f"{100 * res['validity']['te_in_covered_frac']:.0f}"
                          "%, below it the taper extrapolates from the "
                          "measured floor)")
            else:
                route += ", sub-table rates = flat+taper stand-in"
        if getattr(self, "_beam2_shift_s", 0.0) > 0.0:
            route += (f", 2-pulse beam (second at "
                      f"{self.p['beam2_t0_ns']:g} ns)")
        if self._mw_on:
            v = res["validity"]
            t0, t1 = v["mw_window_ns"]
            route += (f", + MW window {v['mw_power_W']:g} W inj "
                      f"{t0:g}-{t1:g} ns (E/N cap {v['mw_en_max_td']:g} Td"
                      f", swarm <eps> at cap "
                      f"{v['mw_eps_swarm_at_cap_eV']:.2f} eV, nu/omega "
                      f"{v['mw_nu_over_omega_at_cap']:.0f}"
                      + ("" if v["mw_cap_in_table"] else
                         "; cap OUTSIDE the swarm table") + ")")
        worst = max((abs(c["residual_rel"]) for c in res["closure"][1:]),
                    default=0.0)
        return (f"mode={self.mode} drive={self.p['drive']} grid="
                f"{self.shape[0]}x{self.shape[1]}x{self.shape[2]} "
                f"({self.n_cells} cells, {len(self.net)} reactions) "
                f"{self.p['pressure_atm']:g} atm | {route} | "
                f"drive {res['final_channels_J']['drive_in']:.3e} J, laser "
                f"{res['final_channels_J']['laser_out']:.3e} J "
                f"(eta {100 * res['eta_intrinsic']:.2f}%) | ledger residual "
                f"max {worst:.2e} | validity: ionization "
                f"{res['validity']['ionization_frac_max']:.1e} of N, F2 "
                f"burn {100 * res['validity']['f2_burn_frac']:.1f}%, "
                f"<eps> in-table "
                f"{100 * res['validity']['te_in_table_frac']:.0f}% | "
                f"wall {res['wall_s']:.1f} s | "
                "cells are independent reactors (no inter-cell transport; "
                "coarse device grid, not a PIC grid)")


def _has_lxcat(gas: str) -> bool:
    from .cross_sections import GAS_FILES, gas_file
    if gas not in GAS_FILES:
        return False
    try:
        gas_file(gas)
        return True
    except FileNotFoundError:
        return False
