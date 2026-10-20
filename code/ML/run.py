"""ML command line: fit LIGHT's channels to the measurements.

    python3 ML/run.py data          # the dataset registry: role, points, source
    python3 ML/run.py consistency   # which datasets agree (leave-one-out)
    python3 ML/run.py fit           # least squares + MCMC on the train set
    python3 ML/run.py nn            # the neural-network estimator, tested
    python3 ML/run.py score         # every set against every dataset
    python3 ML/run.py electron      # cross sections vs swarm data, the EEDF,
                                    # and what they do at Tahoma (about 5 min)
    python3 ML/run.py formation     # the KrF formation chain against Tahoma's
                                    # g0, Isat and g0*Isat/Pin (about 2 min)
    python3 ML/run.py pure          # the F2-free half of the chain against
                                    # pure Kr, Ar and Ne data (about 1 min)
    python3 ML/run.py all           # all of the above (about 8 min)
    python3 ML/run.py tail          # regenerate LIGHT's two-temperature tail
                                    # table (about 5 min; not part of 'all':
                                    # it rewrites a pinned LIGHT data file)
LIGHT_Kr, pure krypton (ml/kr_*.py; LIGHT_KR_VERSION, default v0_3):
    python3 ML/run.py kr-data         # krypton's swarm datasets and roles
    python3 ML/run.py kr-consistency  # which agree (leave-one-out, pruned)
    python3 ML/run.py kr-fit          # least squares, MCMC, neural estimator
    python3 ML/run.py kr-train        # the fit, and the factors file LIGHT_Kr
                                      # pins (not part of 'all': it writes a
                                      # LIGHT_Kr data file)
The first kr command simulates the true solver for the emulator (about
30 min on 6 workers, KR_WORKERS) and caches it under ML/results/.
Results also go to ML/results/<command>.txt (gitignored).
"""
from __future__ import annotations

import argparse
import contextlib
import io
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ml import data, fit, forward, nn, params, score  # noqa: E402


def cmd_data():
    for role in data.ROLES:
        print(f"\n{role.upper()}")
        for d in data.by_role(role):
            err = (f"{d.sigma * 100:.0f} %" if d.sigma_is_fraction
                   else f"+-{d.sigma:.3g}")
            print(f"  {d.name:18s} {d.kind:16s} {len(d.points):3d} points,"
                  f" error {err}")
            print(f"      {d.source}")
            if d.note:
                print(f"      note: {d.note}")


def _groups():
    g = {"Eden": data.by_name("eden_kr", "eden_f2"),
         "Kansas State": data.by_name("ksu_fig1", "ksu_fig2"),
         "Hakuta": data.by_name("hakuta_fig2", "hakuta_fig3"),
         "Yu & Setser": data.by_name("ys_ir300")}
    return g


def cmd_consistency():
    g = _groups()
    print("Fits with each group left out; chi2 per point for every group "
          "(about 1 = agrees with its stated error):")
    head = "  " + f"{'left out':14s}" + "".join(f"{k:>14s}" for k in g) \
        + "      k_BC   K_eq  tau_C  B share"
    print(head)
    for row in fit.consistency(g):
        c = row["consts"]
        print("  " + f"{str(row['left_out'] or 'none'):14s}"
              + "".join(f"{row['chi2'][k]:14.2f}" for k in g)
              + f"   {c.k_BC:.1e}  {c.K_eq:5.2f}  {c.tau_C * 1e9:4.0f}  "
              f"{forward.b_share_1atm(c):6.2f}")


def _print_summary(s):
    for n in params.NAMES:
        lo, mid, hi = s[n]
        if n.startswith("tau"):
            print(f"  {n:8s} {mid * 1e9:7.2f} ns  [{lo * 1e9:.2f}, "
                  f"{hi * 1e9:.2f}]")
        elif n == "K_eq":
            print(f"  {n:8s} {mid:7.2f}     [{lo:.2f}, {hi:.2f}]")
        else:
            print(f"  {n:8s} {mid:9.2e}  [{lo:.2e}, {hi:.2e}]")
    lo, mid, hi = s["b_share_1atm"]
    print(f"  lasing B share at 1 atm, 300 K: {mid:.2f} [{lo:.2f}, {hi:.2f}]"
          " (LIGHT now "
          f"{forward.b_share_1atm(params.light_now()):.2f})")


def cmd_fit(n_steps=2500, device=True):
    train = data.by_role("train")
    print("train set: " + ", ".join(d.name for d in train))
    x = fit.fit(train)
    chain, acc = fit.mcmc(train, x, n_steps=n_steps)
    print(f"MCMC, 40 walkers x {n_steps} steps (acceptance {acc:.2f}): "
          "median [16th, 84th percentile]")
    _print_summary(fit.summarize(chain))
    xm = __import__("numpy").median(chain, axis=0)
    chi = fit.chi2_per_point(xm, [d for d in data.all_datasets()
                                  if d.kind not in data.DEVICE_KINDS])
    print("  chi2 per point at the median: " + "  ".join(
        f"{k} {v:.2f}" for k, v in chi.items()))
    if device:
        lo, mid, hi = score.device_posterior(chain)
        print(f"  Tahoma g0 over 40 posterior draws: {mid:.2f} %/cm "
              f"[{lo:.2f}, {hi:.2f}], measured 6.39 (a score, never a fit "
              "target)")
    return chain


def cmd_nn():
    rng_draws = __import__("numpy").random.default_rng(11)
    np_ = __import__("numpy")
    for label, sets in (("Network 1, train set", data.by_role("train")),
                        ("Network 2, train set + Hakuta (conflicting)",
                         data.by_role("train")
                         + data.by_name("hakuta_fig2", "hakuta_fig3"))):
        t0 = time.time()
        model = nn.train(sets)
        m, s = nn.apply(model)
        print(f"\n{label} ({time.time() - t0:.0f} s)")
        print("  channel    held-out R^2  coverage   answer on the real data")
        for i, n in enumerate(params.NAMES):
            unit, u = (1e9, " ns") if n.startswith("tau") else (1.0, "")
            print(f"  {n:8s}   {model['r2'][i]:10.2f}  {model['coverage'][i]:8.2f}"
                  f"   {math.exp(m[i]) * unit:.3g}{u} [{math.exp(m[i] - s[i]) * unit:.3g}, "
                  f"{math.exp(m[i] + s[i]) * unit:.3g}]")
        draws = np_.clip(m + s * rng_draws.standard_normal((2000, len(m))),
                         nn.P_LO, nn.P_HI)
        shares = [forward.b_share_1atm(params.consts(x)) for x in draws]
        lo, mid, hi = np_.percentile(shares, [16, 50, 84])
        print(f"  lasing B share at 1 atm: {mid:.2f} [{lo:.2f}, {hi:.2f}]")
        chi = fit.chi2_per_point(np_.clip(m, params.LO, params.HI),
                                 data.by_role("train") + data.by_name(
                                     "hakuta_fig2", "hakuta_fig3"))
        print("  chi2 per point of its answer: " + "  ".join(
            f"{k} {v:.2f}" for k, v in chi.items()))


def cmd_score():
    sets = score.candidate_sets()
    rows, dev = score.scorecard(sets)
    names = list(sets)
    print("Constant sets (k_e re-derived for each):")
    for n, c in sets.items():
        print(f"  {n:16s} k_BC {c.k_BC:.1e}  K_eq {c.K_eq:.2f}  tau_C "
              f"{c.tau_C * 1e9:.0f} ns  tau_B {c.tau_B * 1e9:.2f} ns  k_e "
              f"{c.k_e:.1e}  B share 1 atm {forward.b_share_1atm(c):.2f}")
    print("\n  " + f"{'dataset (role)':34s}" + "".join(f"{n:>18s}"
                                                        for n in names))
    for ds, vals in rows:
        if ds.kind in data.DEVICE_KINDS:
            meas = ds.values[0]
            cells = "".join(f"{vals[n]:10.3g} (x{vals[n] / meas:.2f})"
                            for n in names)
        else:
            cells = "".join(f"{vals[n]:17.1f}%" for n in names)
        print(f"  {ds.name + ' (' + ds.role + ')':34s}{cells}")


def cmd_electron():
    from ml import eedf, swarm
    sets = swarm.load()
    by = {s.name: s for s in sets}
    print("1. LoKI-MC cross sections against swarm measurements (two-term "
          "Boltzmann solver;\n   chi2 per point on log errors, about 1 = "
          "agrees within the stated error)")
    print(f"   {'dataset':20s} {'role':9s} {'pts':>4s} {'SIGLO/Biagi':>12s} "
          f"{'TRINITI':>8s}  model/data")
    for ds in sets:
        m = swarm.predict(ds)
        tri = (f"{swarm.chi2(ds, f2_set='F2_TRINITI'):8.2f}"
               if ds.name.startswith("f2") else f"{'':8s}")
        r = m / ds.y
        print(f"   {ds.name:20s} {ds.role:9s} {len(ds.x):4d} "
              f"{swarm.chi2(ds):12.2f} {tri}  {r.min():.2f}-{r.max():.2f}")
    print("\n2. F2 attachment: one scale on its main resonance (0-2.5 eV), "
          "fitted per subset")
    f2 = [n for n in by if n.startswith("f2_katt") and by[n].role == "train"]
    for use in [[n] for n in f2] + [list(swarm.TRAIN_F2), f2]:
        sl, e = swarm.fit_low([by[n] for n in use])
        chi = "  ".join(f"{n[8:]} {swarm.chi2(by[n], sl, 1.0):.2f}"
                        for n in f2)
        mark = "  <- train set" if tuple(use) == swarm.TRAIN_F2 else ""
        print(f"   {'+'.join(n[8:] for n in use):32s} s_low {sl:5.3f} "
              f"(+-{e:.3f} in ln)   chi2: {chi}{mark}")
    s_low, s_err = swarm.fit_low([by[n] for n in swarm.TRAIN_F2])
    c = score.candidate_sets()["fit (train set)"]
    dev = forward.run_device(c)
    b = eedf.beam_eedf(dev)
    b0 = eedf.beam_eedf(dev, ee_on=False)
    print("\n3. The EEDF at Tahoma (LIGHT's 200-400 ns state, fit set). "
          "Rates over the Boltzmann EEDF\n   against the Maxwellian with the "
          "same mean energy:")
    print(f"   beam-sustained with e-e: n_e {b['ne']:.2e} cm^-3, <eps> "
          f"{b['mean']:.2f} eV (LIGHT: n_e {dev['ne']:.2e}, 3/2 Te "
          f"{1.5 * dev['te']:.2f})")
    print("     " + "  ".join(f"{r} {v:.2g}" for r, v in b["f"].items()))
    print(f"   an attached electron carries off {b['e_att']:.2f} eV "
          f"(LIGHT books 3/2 Te = {1.5 * dev['te']:.2f} eV)")
    print(f"   without e-e (a swarm code's physics): n_e {b0['ne']:.2e}, "
          f"attachment x{b0['f']['ATT']:.2f}, G1 x{b0['f']['G1']:.1e}")
    print("\n4. Tahoma (device data: a score, never a fit target). "
          "g0 measured 6.39 %/cm, formation 0.27")
    for name in ("LIGHT now", "fit (train set)"):
        cs = score.candidate_sets()[name]
        rows = [("Maxwellian, (3/2)Te capture (LIGHT today)",
                 forward.run_device(cs))]
        for lab, kw in (("beam EEDF with e-e, exact capture", {}),
                        (f"  + trained attachment x{s_low:.3f}",
                         dict(att=(s_low, 1.0))),
                        ("  bracket: Chen's scale", dict(att=(swarm.fit_low(
                            [by["f2_katt_chen"]])[0], 1.0))),
                        ("  bracket: Trainor & Jacob 500 K scale",
                         dict(att=(swarm.fit_low([by["f2_katt_tj500"]])[0],
                                   1.0)))):
            rows.append((lab, eedf.self_consistent(cs, **kw)[0]))
        print(f"   {name}")
        for lab, d in rows:
            print(f"     {lab:44s} g0 {d['g0']:5.2f} (x{d['g0'] / 6.39:.2f})"
                  f"  formation {d['formation']:.3f} "
                  f"(x{d['formation'] / 0.27:.2f})  n_e {d['ne']:.2e}"
                  f"  Te {d['te']:.2f}")


def cmd_formation():
    from ml import formation, formation_data
    sets = score.candidate_sets()
    print("KrF formation at Tahoma (device data: a score, never a fit "
          "target).\nMeasured: g0 6.39 %/cm, Isat 1.4 +- 0.1 MW/cm2, "
          "g0*Isat/Pin 0.27 (Spectra Technology 1989, p.3-16)\n")
    print("1. The observable. g0*Isat/Pin from the model's own saturation "
          "(248 nm probe flux,\n   Isat fitted at fixed g0 over 0.5-4 "
          "MW/cm2), with LIGHT's formation chain:")
    for n, c in sets.items():
        o = formation.observables(c)
        print(f"   {n:16s} g0 {o['g0']:5.2f}  Isat {o['isat']:5.2f}  "
              f"g0*Isat/Pin {o['eta']:.3f}  (total B+C formation "
              f"{o['formation']:.3f})")
    c = sets["fit (train set)"]
    base = formation.observables(c)
    print("\n2. Where LIGHT's KrF comes from (fit set, KrF per ion pair): "
          + "  ".join(f"{k} {v:.3f}" for k, v in sorted(
              base["by"].items(), key=lambda x: -x[1]) if v > 0))
    ev = base["events"]
    print(f"   Kr* made per pair: beam 0.775, thermal electrons "
          f"{ev['G1'] + ev['C1a'] + ev['C1b'] + ev['C1c']:.3f}, "
          f"Kr2+ + e {ev['R11p']:.3f}")
    print("\n3. Each link at its measured values, one at a time (fit set):")
    print(f"   {'link':8s} {'value':>9s}  {'g0':>5s} {'Isat':>5s} "
          f"{'eta':>6s}   source (status)")
    print(f"   {'LIGHT':8s} {'':9s}  {base['g0']:5.2f} {base['isat']:5.2f} "
          f"{base['eta']:6.3f}")
    for v in formation_data.VALUES:
        if v.link == "k3_tpow":
            continue
        ch = formation.LIGHT_CHAIN.with_(**{v.link: v.value})
        if v.link.startswith("kr2") or v.link.startswith("k_kr2"):
            ch = ch.with_(k3=1.7e-32)        # the Kr2* links need Kr2*
        o = formation.observables(c, ch)
        print(f"   {v.link:8s} {v.value:9.3g}  {o['g0']:5.2f} "
              f"{o['isat']:5.2f} {o['eta']:6.3f}   {v.source.split(',')[0]}"
              f" ({v.status})")
    meas = formation.LIGHT_CHAIN.with_(**formation_data.MEASURED)
    print("\n4. Every measured link at its measured value "
          "(formation_data.MEASURED):")
    for n, cs in sets.items():
        o = formation.observables(cs, meas)
        print(f"   {n:16s} g0 {o['g0']:5.2f}  Isat {o['isat']:5.2f}  "
              f"g0*Isat/Pin {o['eta']:.3f} (x{o['eta'] / 0.27:.2f})")
    om = formation.observables(c, meas)
    print("   The unmeasured links, one at a time from there (fit set):")
    for lab, kw in formation_data.UNMEASURED:
        o = formation.observables(c, meas.with_(**kw))
        print(f"     {lab:52s} g0 {o['g0']:5.2f}  g0*Isat/Pin "
              f"{o['eta']:.3f} ({o['eta'] - om['eta']:+.3f})")


def cmd_pure():
    from ml import pure, pure_data as P
    print("Verification on pure rare gases: the part of the KrF chain without "
          "F2 (sources in\nml/pure_data.py; abstract = read in the abstract "
          "only).\n")
    print("1. The beam's supply: Kr* + Kr+ made per eV deposited, "
          "(1 + N_ex/N_i) / W")
    print(f"   {'gas':3s} {'W: SF':>7s} {'meas':>5s}  {'N_ex/N_i: SF':>12s} "
          f"{'meas':>5s} {'Kannari':>7s}   {'per eV: SF':>10s} {'meas':>7s} "
          f"{'Kannari':>7s}")
    for g in ("Ne", "Ar", "Kr"):
        x = pure.supply(g)
        wm = f"{x['W_meas']:5.1f}" if "W_meas" in x else f"{'-':>5s}"
        nm = f"{x['nex_meas']:5.2f}" if "nex_meas" in x else f"{'-':>5s}"
        pm = (f"{x['per_eV_meas']:7.4f}" if "per_eV_meas" in x
              else f"{'-':>7s}")
        print(f"   {g:3s} {x['W_sf']:7.2f} {wm}  {x['nex_sf']:12.3f} {nm} "
              f"{x['nex_kannari']:7.2f}   {x['per_eV_sf']:10.4f} {pm} "
              f"{x['per_eV_kannari']:7.4f}")
    print("   (measured: W from IAEA-TECDOC-799 / Mano 2021 / PDG, N_ex/N_i "
          "from Saito 2003, alphas)")
    rows, k15, conv = pure.kr_rates()
    print("\n2. Kr ion chemistry: LIGHT against pure-Kr afterglows")
    for te, light, meas in rows:
        print(f"   Kr2+ + e at Te {te:5.3f} eV: LIGHT {light:.2e}, Shiu & "
              f"Biondi {meas:.2e} cm3/s (x{light / meas:.2f})")
    print(f"   Kr+ + 2Kr: LIGHT {k15:.1e}, Smith et al. {conv[1]:.1e}-"
          f"{conv[0]:.1e} cm6/s")
    print("\n3. LIGHT on pure Kr (no F2): VUV yield, light counted to 20 us "
          "after a 450 ns pulse")
    ideal = pure.ideal_vuv("Kr")
    print(f"   measured at low pump (Morozov 2008): 0.42 +- 0.05. One photon "
          f"per excitation and ion with\n   the measured W and N_ex/N_i: "
          f"{ideal:.3f}")
    for pw in (3.3e2, 3.3e3, 3.3e4, 3.3e5):
        r = pure.kr_vuv(pw)
        print(f"   {pw / 1e3:7.2f} kW/cm3: VUV {r['vuv']:.3f}  (LIGHT's "
              f"lossless ceiling {r['ceiling']:.3f}; thermal-electron Kr* "
              f"{r['thermal']:.3f}/pair)")
    r = pure.kr_vuv(3.3e2, thermal=False)
    print(f"   0.33 kW/cm3 without thermal-electron excitation: VUV "
          f"{r['vuv']:.3f}")
    print(f"\n4. Tahoma's pure-Kr run: absorption {P.TAHOMA_PURE_KR.value} "
          f"%/cm at 257.3 nm. Not usable yet:\n   the Kr2+ and Kr2* cross "
          f"sections there are theory or assumption.")


def cmd_tail():
    from ml import tail
    print("The two-temperature tail: r = T_tail / Te from the kinetic EEDF at "
          "LIGHT's own states\n(ml/tail.py). Writes "
          "LIGHT/v0/data/kr_tail_ratio.txt.\n")
    pts = tail.scan()
    path, md5, check = tail.write_table(pts)
    tail.write_states(pts, HERE / "results" / "tail_states.csv")
    print(f"\nwrote {path.name}, md5 {md5}. Pin it in LIGHT/v0/xsec.py "
          "(TAIL_FILE).")
    print("pressure check, table r / kinetic r: "
          + ", ".join(f"{t:.0f} Torr {q:.3f}" for t, q in check))


def _kr_names(rounds):
    return rounds["kept"]


def cmd_kr_data():
    from ml import kr_data
    print("Krypton's electron swarm datasets (ml/kr_data.py): role, kind, "
          "points, E/N range, assumed error\n")
    for d in kr_data.load():
        print(f"  {d.name:20s} {d.role:8s} {d.kind:6s} {len(d.x):3d} pts "
              f"{d.x.min():9.4g}-{d.x.max():<8.4g} Td  {d.sigma:.0%}  "
              f"Tg {d.tg:.0f} K")
        print(f"      {d.source}\n      via {d.via.split(' (')[0]}; "
              f"{d.note}")


def cmd_kr_consistency():
    from ml import kr_train
    s = kr_train.setup()
    c = kr_train.run_consistency(s)
    print("\nConsistency (leave one out; the worst dropped while its chi2 "
          "per point exceeds 3):")
    for i, r in enumerate(c["rounds"], 1):
        print(f"  round {i}: " + ", ".join(
            f"{n} {r['loo'][n]:.2f}" for n in r["kept"]))
    print("\ndropped, in order: " + ", ".join(
        f"{n} ({v:.1f})" for n, v in c["dropped"]))
    print("train set: " + ", ".join(c["kept"]))


def cmd_kr_fit(write=False):
    from ml import kr_light, kr_train
    s = kr_train.setup()
    c = kr_train.run_consistency(s)
    f = kr_train.run_fit(s, c["kept"])
    kr_train.report(s, f)
    if write:
        xs = kr_light.module("xsec")
        path = xs.DATA_DIR / xs.TRAINED_FILE[0]
        rec = kr_train.record(f, "ML/run.py kr-train: least squares on the "
                              "consistency check's train set, the emulator "
                              "checked against the true solver", s, c)
        md5 = kr_swarm_write(rec, path)
        print(f"\nwrote {path}, md5 {md5}: pin it in LIGHT_Kr/"
              f"{kr_light.LIGHT_KR_VERSION}/xsec.py (TRAINED_FILE)")


def kr_swarm_write(rec, path):
    from ml import kr_swarm
    return kr_swarm.write_trained(rec, path)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ML/run.py", description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("command", choices=("data", "consistency", "fit", "nn",
                                        "score", "electron", "formation",
                                        "pure", "tail", "all", "kr-data",
                                        "kr-consistency", "kr-fit",
                                        "kr-train"))
    ap.add_argument("--steps", type=int, default=2500, help="MCMC steps")
    a = ap.parse_args(argv)
    cmds = (["data", "consistency", "fit", "nn", "score", "electron",
             "formation", "pure"] if a.command == "all" else [a.command])
    out_dir = HERE / "results"
    out_dir.mkdir(exist_ok=True)
    for c in cmds:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            {"data": cmd_data, "consistency": cmd_consistency,
             "fit": lambda: cmd_fit(a.steps), "nn": cmd_nn,
             "score": cmd_score, "electron": cmd_electron,
             "formation": cmd_formation, "pure": cmd_pure,
             "tail": cmd_tail, "kr-data": cmd_kr_data,
             "kr-consistency": cmd_kr_consistency, "kr-fit": cmd_kr_fit,
             "kr-train": lambda: cmd_kr_fit(write=True)}[c]()
        text = buf.getvalue()
        print(f"===== {c}\n{text}")
        (out_dir / f"{c}.txt").write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
