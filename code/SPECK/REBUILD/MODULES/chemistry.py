"""RHS assembled from the network — nothing hand-written, everything routed.

This is defect class D2's structural prevention in code: the 872-line
hand-written ``_rhs_u_body.py`` of phantom11 is replaced by a compiler that
walks ``network.py``'s data rows and assembles a vectorized right-hand side
over an arbitrary number of device-grid cells. A term cannot hide because no
term is written by hand; a sign cannot flip silently because stoichiometry
comes from the (balance-checked) row itself.

ENERGY ROUTING (the ledger's feed, all per event, computed at compile time
from the row + the species registry — see network.py docstring for rules):

  u_from_e   drawn from the electron gas (endothermic e-impact + activation)
  u_to_e     added to the electron gas (Penning-class ejection)
  u_to_gas   released to (or drawn from, if negative) the gas
  u_rad      emitted quantum of a radiative row (escapes, except the upper
             laser level's K_axial-captured fraction — split by the engine)

plus the runtime term for electron-consuming rows (attachment/DR): the
captured electron's mean thermal energy (3/2)Te moves electron gas → gas.

The photon-species coupling (stimulated emission on the upper level, the
absorber inventory with PHOTO_PRODUCTS routing) lives here too, because it
transforms species; the global photon equation itself is the engine's
(laser.photon_rhs on gain-volume averages).
"""

from __future__ import annotations

import numpy as np

from .network import PHOTO_PRODUCTS
from .species import charge_of, e_form_of
from .units import C_LIGHT, E_CHARGE, EPS0, EV, K_B


def ion_ion_screening(ne, n_anion, tg):
    """Debye screened-capture factor for ion-ion recombination (0.1.14,
    chain 3 of the amplification audit).

    Kannari 1985 p.4313 (via Morgan [35]): plasma shielding cuts ion-ion
    recombination at pulse ion densities — he prints 0.397 at 1e14 cm^-3
    pairs in 2-atm Ar and applied it as a CONSTANT. Here the
    parameter-free screened-capture form is used instead:

        S = 1 / (1 + R_c/lambda_D),   R_c = 2 e^2 / (3 * 4 pi eps0 kB Tg)

    (R_c: the radius where |U_Coulomb| = (3/2) kT), lambda_D over ALL
    charge carriers (ne + cations + anions = 2 (ne + n_anion), the slow
    species doing the screening) at Tg. At his printed point this form
    gives 0.695 vs Morgan's Monte-Carlo 0.397 — a MILDER cut, x1.75
    apart, RECORDED (Morgan's paper is not in-repo; the form's spread is
    part of the error budget). S -> 1 as the plasma vanishes. Applied to
    the flannery/lee_ii rate kinds by CompiledNetwork.rhs."""
    n_scr = np.maximum(2.0 * (ne + n_anion), 1.0)
    tg_k = np.maximum(tg, 1.0)
    lam_d = np.sqrt(EPS0 * K_B * tg_k / (n_scr * E_CHARGE ** 2))
    r_cap = 2.0 * E_CHARGE ** 2 / (3.0 * 4.0 * np.pi * EPS0 * K_B * tg_k)
    return 1.0 / (1.0 + r_cap / lam_d)

_TE_MIN_EV = 0.026          # rate-evaluation floor (room temperature)
VIB_EV = 0.2                # F2 vibrational quantum per 'vib' event
                            # (REPO: spectre engine.py:258)
# electron/atom mass ratios for the elastic-drag term (2 m_e/M)
_ME_AMU = 5.485799e-4
_GAS_AMU = {"Ne": 20.1797, "Ar": 39.948, "Kr": 83.798, "F2": 37.9968,
            "Xe": 131.293}
_KB_EV_K = 8.617333262e-5


def petrov17_f2_k(channel: str, P, x_f2):
    """Beam-sustained-EEDF rate coefficients for the low-threshold e+F2
    channels — Petrov, Wolford, Petrova, Giuliani, Obenschain, JAP 122,
    133301 (2017), Table I (read first-hand 2026-09-01), the fits to
    their beam-deposition EEDF solutions (their Fig. 6; the memo
    petrova2022 Table A2 k23/k24/k25 carries the same formulas).

    ``P`` is the normalized power deposition 2*P_beam[MW/cm^3]/p[atm]
    + 1e-3 (their stability floor), ``x_f2`` the F2 number fraction.
    Returns k in m^3/s (the paper prints cm^3/s). Applied with the BULK
    electron density: the fits fold the beam-secondary tail's shape into
    the coefficient, replacing (never adding to) the <eps>-indexed swarm
    coefficient for the same channel — one owner per pathway. Domain:
    their computed grid spans P ~ 0.01-1+, x_f2 <= 1%; the +1e-3 floor
    keeps the ramps finite."""
    P = np.maximum(np.asarray(P, dtype=float), 1e-3)
    x = np.maximum(np.asarray(x_f2, dtype=float), 0.0)
    ln_p = np.log(P)
    rx = np.sqrt(x)
    if channel == "vib":
        a = 1.0e-9
        b = (2.0 - 0.45 * ln_p) * 1.0e-9
        c = 50.0 / (1.0 + 4.0 * P)
        k_cm3 = a + b / (1.0 + c * rx)
    elif channel == "att":
        b = (2.3 - 0.5 * ln_p) * 1.0e-9
        c = 50.0 / (1.0 + 4.0 * P)
        k_cm3 = b / (1.0 + c * rx)
    elif channel == "dis":
        a = (0.7 + 1.7 * P ** 0.35) * 1.0e-9
        b = 0.2 + 0.15 * P ** 0.4
        c = 1.0e3 / (1.0 + 4.0 * P)
        k_cm3 = a / (1.0 + 1.0 / (b + c * x) + 10.0 * rx)
    else:
        raise ValueError(f"unknown petrov17 channel {channel!r}")
    return k_cm3 * 1.0e-6


class CompiledNetwork:
    """The network compiled against a StateLayout + electron-rate tier."""

    def __init__(self, rows, layout, tiers, optics, petrov17=False):
        self.rows = list(rows)
        self.lay = layout
        self.tiers = tiers
        self.optics = dict(optics)
        self.upper = optics["upper"]
        self.hv_eV = None       # set below from the cavity wavelength
        i_upper = layout.index.get(self.upper)
        if i_upper is None:
            raise KeyError(f"upper level {self.upper!r} not in layout")
        self._compiled = []
        for r in self.rows:
            heavy = []
            for s in r.heavy_reactants():
                if s not in layout.index:
                    raise KeyError(f"{r.rid}: species {s!r} not in this "
                                   "mode's layout")
                heavy.append(layout.index[s])
            stoich = {}
            for s in r.reactants:
                if s in layout.index:
                    stoich[s] = stoich.get(s, 0) - 1
            for s in r.products:
                if s in layout.index:
                    stoich[s] = stoich.get(s, 0) + 1
            pairs = [(layout.index[s], float(c)) for s, c in stoich.items()
                     if c != 0]
            ne_pow = sum(1 for s in r.reactants if s == "e")
            release = r.release_eV()
            e_cost = float(r.e_cost_override_eV or 0.0)
            u_from_e = u_to_e = u_to_gas = u_rad = 0.0
            is_upper_rad = False
            if r.kind == "rad":
                hv = float(r.hv_eV) if r.hv_eV is not None else max(release, 0.0)
                if r.energy_to == "gas":
                    u_to_gas += release
                else:
                    u_rad += hv
                    u_to_gas += release - hv
                    if r.heavy_reactants()[0] == self.upper:
                        is_upper_rad = True
            elif r.is_electron_row():
                if release < 0.0:
                    u_from_e += -release
                elif r.energy_to == "electrons":
                    u_to_e += release
                else:
                    u_to_gas += release
                u_from_e += e_cost
                u_to_gas += e_cost
            else:
                if release < 0.0:
                    u_to_gas += release
                else:
                    dest = r.energy_to or ("electrons" if r.ne_delta() > 0
                                           else "gas")
                    if dest == "electrons":
                        u_to_e += release
                    elif dest == "radiation":
                        u_rad += release
                    else:
                        u_to_gas += release
            capture = max(-r.ne_delta(), 0)
            is_ii = r.rate[0] in ("flannery", "lee_ii")
            self._compiled.append((r, heavy, ne_pow, pairs, u_from_e, u_to_e,
                                   u_to_gas, u_rad, capture, is_upper_rad,
                                   float(r.tg_pow), is_ii))
        # ion-ion Debye screening needs the anion inventory (cations are
        # ne + anions by quasineutrality; total charge carriers screen)
        self._anion_idx = [i for i, s in enumerate(self.lay.species)
                           if charge_of(s) < 0]
        # --- photon-absorber routing, precompiled --------------------------
        self._absorbers = []
        for sp, sig in self.optics["absorbers"].items():
            i_sp = layout.index.get(sp)
            if i_sp is None:
                continue
            if sp not in PHOTO_PRODUCTS:
                raise KeyError(f"absorber {sp!r} has no PHOTO_PRODUCTS "
                               "routing — an absorbed photon must transform "
                               "its absorber (ledger rule)")
            prods = PHOTO_PRODUCTS[sp]
            heavy_prods = [p for p in prods if p != "e"]
            for p in heavy_prods:
                if p not in layout.index:
                    raise KeyError(f"photo-product {p!r} of {sp!r} not in "
                                   "layout")
            de_pot = (sum(e_form_of(p) for p in prods) - e_form_of(sp))
            q_heavy = (sum(charge_of(p) for p in heavy_prods)
                       - charge_of(sp))
            n_freed = q_heavy          # +1 -> a free electron was ejected
            self._absorbers.append((i_sp, float(sig),
                                    [layout.index[p] for p in heavy_prods],
                                    de_pot, n_freed > 0, sp))
        self._i_upper = i_upper
        self._frag_idx = [layout.index[f] for f in self.optics["fragments"]]
        # --- petrov17 beam-EEDF e+F2 rows ([0.1.36], flag-gated) -----------
        # Marked by PATTERN (reactants e+F2, loki-rated, the three
        # low-threshold buckets), immune to rid drift and the repulsive
        # reroute. The vib sink is handled at its own block in rhs.
        self._petrov17 = bool(petrov17)
        self._i_f2 = layout.index.get("F2")
        self._p17_att = set()
        self._p17_dis = {}                     # compiled idx -> bucket
        if self._petrov17 and self._i_f2 is not None:
            for i, entry in enumerate(self._compiled):
                r = entry[0]
                if (set(r.reactants) == {"e", "F2"}
                        and r.rate[0] == "loki"):
                    if r.rate[1] == "F-":
                        self._p17_att.add(i)
                    elif r.rate[1] in ("F2(a3)", "F2(A1)"):
                        self._p17_dis[i] = r.rate[1]

    def set_photon_energy(self, hv_eV: float):
        self.hv_eV = float(hv_eV)
        up = e_form_of(self.upper)
        frag = sum(e_form_of(self.lay.species[i]) for i in self._frag_idx)
        self._stim_gas_eV = up - frag - hv_eV
        if self._stim_gas_eV < -1e-9:
            raise ValueError(
                f"stimulated emission would create energy: E_form({self.upper})"
                f" = {up} eV < fragments {frag} + hv {hv_eV}")
        # The upper level's SPONTANEOUS row must radiate the EXACT cavity
        # quantum (its hand-stated hv_eV is 5.00 vs the exact 4.9997 at
        # 248 nm) or the ledger leaks 0.3 meV/event — caught by the
        # per-event identity gate 2026-08-04. Recompute its gas remainder
        # from the exact quantum here.
        for i, entry in enumerate(self._compiled):
            (r, heavy, ne_pow, pairs, u_from_e, u_to_e, u_to_gas, u_rad,
             capture, is_upper_rad, tg_pow, is_ii) = entry
            if is_upper_rad:
                u_to_gas = r.release_eV() - self.hv_eV
                if u_to_gas < -1e-9:
                    raise ValueError(
                        f"{r.rid}: exact quantum {self.hv_eV} eV exceeds "
                        f"the row's release {r.release_eV()} eV")
                self._compiled[i] = (r, heavy, ne_pow, pairs, u_from_e,
                                     u_to_e, u_to_gas, u_rad, capture,
                                     is_upper_rad, tg_pow, is_ii)
        for (_i, _s, _p, de_pot, freed, sp) in self._absorbers:
            if hv_eV - de_pot < -1e-9:
                raise ValueError(f"absorber {sp}: potential change {de_pot} eV"
                                 f" exceeds the photon {hv_eV} eV")

    # ------------------------------------------------------------------ rates
    def _tier_k(self, bucket, te, en_td):
        """LOKI/KINETIC lookup: by E/N when the circuit drive supplies one
        (the tables' NATIVE index — no <eps> mapping), else by Te."""
        if en_td is not None:
            return self.tiers.k_by_en(bucket, float(en_td))
        return self.tiers.k(bucket, np.maximum(te, _TE_MIN_EV))

    def _k(self, spec, te, ntot, en_td=None):
        kind = spec[0]
        if kind == "const":
            return spec[1]
        if kind == "te_fit":
            A, B, C = spec[1], spec[2], spec[3]
            t = np.maximum(te, _TE_MIN_EV)
            return A * t ** B * np.exp(C / t)
        if kind == "stepwise":
            t = np.maximum(te, _TE_MIN_EV)
            A = spec[2] if len(spec) > 2 else 1.56e-13
            return A * t ** 0.71 * np.exp(-spec[1] / t)
        if kind == "dr_sqrt":
            return spec[1] / np.sqrt(np.maximum(te, 0.02))
        if kind == "loki":
            return self._tier_k(spec[1], te, en_td)
        if kind == "flannery":
            n_am = ntot / _amagat()
            return (spec[1] * n_am / (n_am * n_am + 1.7)) * spec[2] * 1e-6
        if kind == "lee_ii":
            n_am = ntot / _amagat()
            A, b, c = spec[1], spec[2], spec[3]
            return (A * n_am / (n_am * n_am + b * n_am + c)) * 1e-6
        raise ValueError(f"unknown rate kind {kind!r}")

    def _p17_ctx(self, dens, te, ntot, en_td, p17_norm):
        """Per-eval petrov17 override context, or None when inactive.

        Returns (x_f2 (C,), k_att (C,), k_dis_total (C,), w bucket->(C,),
        p17_norm): the beam-EEDF coefficients for this state, with the
        total neutral-dissociation coefficient branched across the marked
        repulsive rows IN PROPORTION to their own swarm rates at this Te
        (the measured branching shape is kept; the magnitude is the
        beam-EEDF total; equal split where the swarm rates vanish)."""
        if p17_norm is None or not self._petrov17 or self._i_f2 is None:
            return None
        x = dens[self._i_f2] / max(float(ntot), 1.0)
        k_att = petrov17_f2_k("att", p17_norm, x)
        k_dis = petrov17_f2_k("dis", p17_norm, x)
        w = {}
        if self._p17_dis:
            buckets = sorted(set(self._p17_dis.values()))
            kb = {}
            for b in buckets:
                v = np.asarray(self._tier_k(b, te, en_td), dtype=float)
                kb[b] = np.maximum(np.broadcast_to(v, x.shape), 0.0)
            tot = np.zeros(x.shape)
            for b in buckets:
                tot = tot + kb[b]
            eq = 1.0 / len(buckets)
            for b in buckets:
                w[b] = np.where(tot > 0.0, kb[b] / np.maximum(tot, 1e-300),
                                eq)
        return (x, k_att, k_dis, w, p17_norm)

    # -------------------------------------------------------------------- rhs
    def row_fluxes(self, dens, te, tg, ne, ntot, en_td=None,
                   p17_norm=None) -> dict:
        """Per-row reaction flux [m^-3 s^-1] (rid -> (C,) array) at one
        state — the SAME k x heavy densities x n_e^pow x T_g / screening
        arithmetic the rhs loop applies, exposed for diagnostics (channel
        budgets: which rows form / destroy a species, integrated over a
        run). Reaction rows only: the stimulated-emission, ASE, absorber
        and elastic/vib terms are not rows and are not included. Gated:
        sum(coeff x row flux) == rhs ddens at n_ph = 0, ase_rate = 0.
        """
        S, C = dens.shape
        tg_fac = (np.maximum(tg, 1.0) / 300.0)
        n_anion = (dens[self._anion_idx].sum(axis=0)
                   if self._anion_idx else np.zeros(C))
        s_scr = ion_ion_screening(ne, n_anion, tg)
        p17 = self._p17_ctx(dens, te, ntot, en_td, p17_norm)
        out = {}
        for i_row, (r, heavy, ne_pow, pairs, u_from_e, u_to_e, u_to_gas,
                    u_rad, capture, is_upper_rad, tg_pow, is_ii) in \
                enumerate(self._compiled):
            k = self._k(r.rate, te, ntot, en_td)
            if p17 is not None:
                if i_row in self._p17_att:
                    k = p17[1]
                elif i_row in self._p17_dis:
                    k = p17[2] * p17[3][self._p17_dis[i_row]]
            if tg_pow != 0.0:
                fl = r.tg_floor          # 0.1.40: a + b Tg^n, normalised
                k = k * (tg_fac ** tg_pow if fl == 0.0
                         else fl + (1.0 - fl) * tg_fac ** tg_pow)
            if is_ii:
                k = k * s_scr
            flux = k if not np.isscalar(k) else np.full(C, float(k))
            flux = np.array(flux, dtype=float, copy=True)
            for i_h in heavy:
                flux *= dens[i_h]
            if r.kind != "rad" and ne_pow:
                flux *= ne ** ne_pow
            np.maximum(flux, 0.0, out=flux)
            out[r.rid] = flux
        return out

    def row_stoichiometry(self) -> dict:
        """rid -> list of (species index, net coefficient) the row applies."""
        return {r.rid: list(pairs) for (r, _h, _n, pairs, *_rest)
                in self._compiled}

    def rhs(self, dens, te, tg, ne, n_ph, ntot, en_td=None,
            ase_rate=0.0, p17_norm=None):
        """All chemistry + photon-species terms for one evaluation.

        Inputs: dens (S, C), te/tg/ne (C,), n_ph scalar, ntot scalar.
        Returns dict with ddens (S, C) [m^-3 s^-1], du_e / dgas / rad_escape
        [W/m^3], and the per-cell photon rates r_ste / r_x / r_se
        [m^-3 s^-1] for the engine's photon equation and ledger split.
        """
        S, C = dens.shape
        ddens = np.zeros((S, C))
        du_e = np.zeros(C)          # W/m^3 into the electron gas
        dgas = np.zeros(C)          # W/m^3 into the gas
        rad = np.zeros(C)           # W/m^3 escaping as spont/VUV radiation
        r_se = np.zeros(C)
        # swarm-residual bookkeeping (E/N-direct path only; the sink block
        # after the elastic/vib section explains and consumes these)
        p_loki_drain = np.zeros(C)  # W/m^3 the table-rate rows drain
        birth_loki = np.zeros(C)    # electrons born /m^3/s by table rows
        tg_fac = (np.maximum(tg, 1.0) / 300.0)
        # ion-ion Debye screened-capture factor (0.1.14, chain 3) — the
        # physics and the Kannari/Morgan corroboration live on
        # ion_ion_screening's docstring
        n_anion = (dens[self._anion_idx].sum(axis=0)
                   if self._anion_idx else np.zeros(C))
        s_scr = ion_ion_screening(ne, n_anion, tg)
        p17 = self._p17_ctx(dens, te, ntot, en_td, p17_norm)
        for i_row, (r, heavy, ne_pow, pairs, u_from_e, u_to_e, u_to_gas,
                    u_rad, capture, is_upper_rad, tg_pow, is_ii) in \
                enumerate(self._compiled):
            k = self._k(r.rate, te, ntot, en_td)
            if p17 is not None:
                if i_row in self._p17_att:
                    k = p17[1]
                elif i_row in self._p17_dis:
                    k = p17[2] * p17[3][self._p17_dis[i_row]]
            if tg_pow != 0.0:
                fl = r.tg_floor          # 0.1.40: a + b Tg^n, normalised
                k = k * (tg_fac ** tg_pow if fl == 0.0
                         else fl + (1.0 - fl) * tg_fac ** tg_pow)
            if is_ii:
                k = k * s_scr
            flux = k if not np.isscalar(k) else np.full(C, float(k))
            flux = np.array(flux, dtype=float, copy=True)
            for i_h in heavy:
                flux *= dens[i_h]
            if r.kind != "rad" and ne_pow:
                flux *= ne ** ne_pow
            if not np.any(flux > 0.0):
                continue
            np.maximum(flux, 0.0, out=flux)
            for i_s, c in pairs:
                ddens[i_s] += c * flux
            f_ev = flux * EV
            if u_from_e or u_to_e:
                du_e += (u_to_e - u_from_e) * f_ev
            if en_td is not None and ne_pow and r.rate[0] == "loki":
                p_loki_drain += (u_from_e - u_to_e) * f_ev
                nb = r.ne_delta()
                if nb > 0:
                    birth_loki += nb * flux
            if u_to_gas:
                dgas += u_to_gas * f_ev
            if capture:
                cap = 1.5 * np.maximum(te, _TE_MIN_EV) * capture * f_ev
                du_e -= cap
                dgas += cap
            if u_rad:
                if is_upper_rad:
                    r_se += flux          # engine splits K-capture vs escape
                else:
                    rad += u_rad * f_ev
        # --- stimulated emission on the upper level ------------------------
        sig_se = self.optics["sigma_se_m2"]
        r_ste = max(n_ph, 0.0) * C_LIGHT * sig_se * dens[self._i_upper]
        ddens[self._i_upper] -= r_ste
        for i_f in self._frag_idx:
            ddens[i_f] += r_ste
        dgas += self._stim_gas_eV * r_ste * EV
        # --- in-loop ASE sink (S6 closure, 2026-08-04): the engine prices
        # a per-upper-molecule stimulated rate from the 4-pi ASE flux of
        # the CURRENT gain field; routing is identical to stimulated
        # emission except the photons ESCAPE (radiated channel) instead of
        # entering the cavity mode -----------------------------------------
        if ase_rate > 0.0:
            r_ase = ase_rate * dens[self._i_upper]
            ddens[self._i_upper] -= r_ase
            for i_f in self._frag_idx:
                ddens[i_f] += r_ase
            rad += self.hv_eV * r_ase * EV
            dgas += self._stim_gas_eV * r_ase * EV
        # --- absorber inventory (photo-product routing) --------------------
        r_x = np.zeros(C)
        for (i_a, sig, prod_idx, de_pot, frees_e, _sp) in self._absorbers:
            r_a = max(n_ph, 0.0) * C_LIGHT * sig * dens[i_a]
            r_x += r_a
            ddens[i_a] -= r_a
            for i_p in prod_idx:
                ddens[i_p] += r_a
            surplus = (self.hv_eV - de_pot) * r_a * EV
            if frees_e:
                du_e += surplus
            else:
                dgas += surplus
        # --- elastic drag + vibrational sink (electron <-> gas) ------------
        tg_ev = _KB_EV_K * tg
        p_el = np.zeros(C)
        for g in ("Ne", "Ar", "Kr", "F2", "Xe"):
            i_g = self.lay.index.get(g)
            if i_g is None:
                continue
            kel = self._tier_k(f"el:{g}", te, en_td)
            if np.isscalar(kel):
                kel = np.full(C, float(kel))
            p_el += ((2.0 * _ME_AMU / _GAS_AMU[g]) * kel * dens[i_g]
                     * 1.5 * (np.maximum(te, _TE_MIN_EV) - tg_ev))
        p_el *= ne * EV                       # W/m^3, e -> gas (signed)
        i_f2 = self.lay.index.get("F2")
        p_vib = np.zeros(C)
        if i_f2 is not None:
            if p17 is not None:
                # beam-EEDF vibrational coefficient (petrov17, same
                # 0.2 eV/event energy routing as the swarm path)
                kvib = petrov17_f2_k("vib", p17[4], p17[0])
            else:
                kvib = self._tier_k("vib", te, en_td)
            if np.isscalar(kvib):
                kvib = np.full(C, float(kvib))
            p_vib = kvib * dens[i_f2] * ne * VIB_EV * EV
        du_e -= p_el + p_vib
        dgas += p_el + p_vib
        # --- swarm-residual sink (E/N-direct path only; 0.1.27) ------------
        # At swarm equilibrium the tables' own power balance reads
        #   e v_d(E/N) E = elastic + ALL inelastic thresholds + growth,
        # but the bucketed rows book only the level-lumped share of the
        # sigma-set's inelastic spend (measured: 12% of the Joule input
        # unbooked at 14.1 Td in the Chen fill, rising to 45% at 150 Td).
        # On this path the bucket rates are E/N-indexed — Te-FLAT — so the
        # unbooked share has no Te-responsive absorber (Ne/Ar elastic
        # alone would need Te ~ 34 eV to soak 12%) and the fluid Te
        # floated x1.67 above the table's own <eps> (the 0.1.25 anchor's
        # x1.55 Te row, attributed 2026-08-26). The residual is DERIVED
        # from the same table (v_d, <eps>) plus this network's own row
        # inventory — nothing fitted: with it, dTe/dt closes EXACTLY at
        # Te = (2/3)<eps>(E/N), with ionization growth credited as
        # dilution (the fluid's built-in mechanism; hence births, not
        # net, in the growth term — attachment's 1.5 Te capture drain
        # cancels in the Te equation). Routed e -> gas like the vib sink
        # (untracked internal states; the gas store keeps the ledger
        # closed). Clamped >= 0: in-table the sigma-set never implies
        # negative unbooked power; the below-floor stand-in region can.
        if (en_td is not None and hasattr(self.tiers, "drift_velocity")
                and hasattr(self.tiers, "mean_energy_at_en")):
            eps_sw = float(self.tiers.mean_energy_at_en(en_td))
            e_vm = en_td * 1e-21 * ntot
            p_in = ne * EV * float(self.tiers.drift_velocity(en_td)) * e_vm
            te_sw = np.full(C, (2.0 / 3.0) * eps_sw)
            p_el_sw = np.zeros(C)
            for g in ("Ne", "Ar", "Kr", "F2", "Xe"):
                i_g = self.lay.index.get(g)
                if i_g is None:
                    continue
                kel = self._tier_k(f"el:{g}", te_sw, en_td)
                if np.isscalar(kel):
                    kel = np.full(C, float(kel))
                p_el_sw += ((2.0 * _ME_AMU / _GAS_AMU[g]) * kel * dens[i_g]
                            * 1.5 * (te_sw - tg_ev))
            p_el_sw *= ne * EV
            p_res = np.maximum(
                p_in - p_loki_drain - p_el_sw - p_vib
                - eps_sw * EV * birth_loki, 0.0)
            du_e -= p_res
            dgas += p_res
        return {"ddens": ddens, "du_e": du_e, "dgas": dgas,
                "rad_escape": rad, "r_ste": r_ste, "r_x": r_x, "r_se": r_se}


def _amagat():
    from .units import AMAGAT_M3
    return AMAGAT_M3
