"""The KrF electron-quench measurement re-read under LIGHT's two-level
KrF(B,C) model (evidence script, not LIGHT code).

What was measured. Trainor & Jacob (Duzy et al., AVCO final report 1980,
OSTI 5268869, pp. 19-23, sources/Duzy_AVCO_RGHlasers_1980_OSTI5268869.pdf;
journal version APL 37, 675 (1980)) recorded the steady B-X sidelight of
e-beam-pumped Kr/F2 (190 and 95 Torr Kr, F2 0.5-4 Torr) against [F2], with
the electron density n_e = S/(k_att [F2]). After subtracting the F2 quench
term, slope over intercept gives the product
    k_e * tau / k_att = 2.9e-7 s                          (p. 21)
They took tau = 6.5 ns and k_att = 4.5e-9 cm3/s, so k_e = 2.0e-7 (p. 21,
Table III, four runs 1.7-2.2e-7).

The value LIGHT uses. Kushner, Hanson & Schneider, APL 55, 2482 (1989)
(docs/Recent_KrF_Papers/kushner_krfs.pdf) re-analysed the same data for
KrF(B) alone:
  - k_a (1-2)e-9 from their EEDF for these conditions (T_e,eff 1.58 eV)
  - Kr2+ recombination added as an electron loss
  - tau = 5.5 ns, which is 6.5 ns shortened by Kr quench and the trimer
They got 6.5e-8 (k_a 2e-9) and 3.7e-8 (k_a 1e-9) and recommend 3-6e-8.
LIGHT's 4.5e-8 is that recommendation, applied to both B and C.

Two-level reading. With B and C coupled, the lifetime in the product is
the Stern-Volmer lifetime of B-X emission (two_level.sv_tau). It depends on
the B/C constants and on whether KrF is formed in B or in C. This script
re-reads the product with that lifetime and with LIGHT's own F2 attachment
rate (vendored SIGLO cross section, Maxwellian), keeping KHS's
recombination factor. Their printed results back out to 0.62 at k_a 2e-9
and 0.70 at 1e-9.

    /usr/bin/python3 LIGHT/docs/krf_data/electron_quench.py     # < 1 s
"""
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
# marker walk to LIGHT/ (the directory holding run_v0.py), no parent counting
_LIGHT = next(p for p in _HERE.parents if (p / "run_v0.py").is_file())
sys.path.insert(0, str(_LIGHT))
from v0 import two_level as L
from v0 import xsec

PRODUCT = 2.9e-7          # s, k_e tau / k_att (TJ p. 21)
KR_TORR = 190.0           # Fig. 9 case, 7.6 A/cm2 (TJ p. 21)
F_REC_KHS = 0.62          # KHS 6.5e-8 / (2.9e-7 x 2e-9 / 5.5e-9)
TE_KHS = 1.58             # eV, KHS's computed T_e,eff for these conditions


def k_e(tau_s, k_att, f_rec=1.0):
    return PRODUCT * k_att * f_rec / tau_s


def sv_tau(c, source=(1.0, 0.0), kr_torr=KR_TORR):
    return L.sv_tau(c, source_B=source[0], source_C=source[1],
                    kr=kr_torr * L.TORR_300K)


def main():
    c0 = L.Consts.from_light()
    katt = xsec.rate_table("F2", "attachment")
    ka_light = float(katt(TE_KHS))
    print(f"TJ product k_e tau / k_att = {PRODUCT:.1e} s at {KR_TORR:.0f} "
          f"Torr Kr. LIGHT's k_att: {ka_light:.2e} at T_e {TE_KHS} eV, "
          f"{float(katt(0.667)):.2e} at 0.667 eV (a 1-eV mean)\n")
    print(f"{'reading':62s} {'tau ns':>6s} {'k_att':>8s} {'f_rec':>5s} "
          f"{'k_e':>8s}")

    def line(label, tau, ka, fr):
        print(f"{label:62s} {tau * 1e9:6.1f} {ka:8.2e} {fr:5.2f} "
              f"{k_e(tau, ka, fr):8.2e}")
    line("TJ as printed (one state, 6.5 ns)", 6.5e-9, 4.5e-9, 1.0)
    line("KHS, k_a 2e-9 (one state, 5.5 ns)", 5.5e-9, 2.0e-9, F_REC_KHS)
    line("KHS, k_a 1e-9", 5.5e-9, 1.0e-9, 0.70)
    tau_b = sv_tau(c0)
    tau_bc = sv_tau(c0, (0.5, 0.5))
    line("LIGHT two-level, KrF formed in B", tau_b, ka_light, F_REC_KHS)
    line("LIGHT two-level, KrF formed half in B, half in C", tau_bc,
         ka_light, F_REC_KHS)
    print("\nHow the two-level value moves (KrF formed in B, other inputs "
          "as the line above):")
    for kbc in (1e-11, 3e-11, 1.2e-10, 3e-10):
        tau = sv_tau(c0.with_(k_BC=kbc))
        line(f"  k_BC {kbc:.1e} (LIGHT {c0.k_BC:.2e})", tau, ka_light,
             F_REC_KHS)
    for te in (0.667, 1.0, 2.0):
        line(f"  T_e {te} eV for k_att (recombination factor kept)", tau_b,
             float(katt(te)), F_REC_KHS)
    line("  no recombination correction", tau_b, ka_light, 1.0)


if __name__ == "__main__":
    main()
