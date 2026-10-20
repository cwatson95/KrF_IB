#!/usr/bin/env python3
"""gases_twin.py -- the gases of lxcat.html as LoKI-MC's own setups define them (Input/*/..._setup.in), for derive.py.

Populations (absolute: a state's share of its gas's molecules):
  * the target of the Elastic / Effective block: 1 (every atom, F2(V0), Cl2, the ramp gas);
  * N2: N2(X,v=*) Boltzmann at Tg with the energies of Nitrogen/N2energyVibLevels.txt (weights 1), N2(X,v=0,J=*)
    Boltzmann at Tg inside v = 0 with E_J = B J(J+1) and LoKI-MC's rotationalDegeneracy_N2, 3(1 + (1 + (-1)^J)/2)(2J+1);
  * O2: O2(X,v=*) Boltzmann with morseOscillatorEnergy (weights 3, cancelling), O2(X,v=0,J=*) with 2J+1;
  * CO: CO(X,v=*) harmonicOscillatorEnergy (weights 1), CO(X,v=0,J=*) with 2J+1;
  * CO2: the setup's explicit v=000 0.9254, v=010 0.0746 (and its statistical weights);
  * H2O: H2O(X,v=*) Boltzmann with H2O/H2O_vibEnergy_Itikawa.txt (the setup's rotational set is not carried);
  * every other state: 0.
A Boltzmann set is normalised over the states the files name (StatePropertyFunctions boltzmannPopulation), and a
child's population is its parent's times its share.  An Effective set's Elastic is EedfGas::elasticFromEffective
CrossSection's: the populations at 300 K (not Tg), or the CO2 file's, subtract every collision (and de-excitation)
on the effective table's energies, and negatives are clipped to zero.
"""
import math, os
import numpy as np
import lx_lib as L

HBAR_EV = 6.626070040e-34 / (2 * math.pi * 1.6021766208e-19)   # Constant::planckReducedInEV

SETS = {
    'Ar': (['Argon/Ar_LXCat.txt'], 'Ar'), 'Kr': (['Krypton/Kr_LXCat.txt'], 'Kr'), 'Xe': (['Xenon/Xe_LXCat.txt'], 'Xe'),
    'He': (['Helium/He_LXCat.txt'], 'He'), 'Ne': (['Neon/Ne_LXCat.txt'], 'Ne'), 'NeB': (['Neon/Ne_LXCat_mt.txt'], 'Ne'),
    'N2': (['Nitrogen/N2_LXCat.txt', 'Nitrogen/N2_rot_LXCat.txt'], 'N2'),
    'O2': (['Oxygen/O2_LXCat.txt', 'Oxygen/O2_rot_LXCat.txt'], 'O2'), 'O': (['Oxygen/O_LXCat.txt'], 'O'),
    'F2': (['Fluorine2/F2_LXCat.txt'], 'F2'), 'Cl2': (['Chlorine/Cl_LXCat.txt'], 'Cl2'),
    'CO': (['CO/iso/CO_LXCat.txt', 'CO/iso/CO_rot_LXCat.txt'], 'CO'), 'CO2': (['CO2/CO2_LXCat.txt'], 'CO2'),
    'H2O': (['H2O/iso/H2O_LXCat.txt'], 'H2O'), 'ramp': (['ReidRampGas/ReidRampGas_LXCat.txt'], None),
}


def state_parts(name):
    """'N2(X,v=0,J=2)' -> ('N2', 'X', '0', '2'); 'Ar(1S0)' -> ('Ar', '1S0', None, None)"""
    i = name.index('(')
    gas, inner = name[:i], name[i + 1:name.rindex(')')]
    f = inner.split(',')
    ele = f[0]
    vib = rot = None
    for t in f[1:]:
        if t.startswith('v='):
            vib = t[2:]
        elif t.startswith('J='):
            rot = t[2:]
    return gas, ele, vib, rot


def boltz(levels, T):
    """levels: {key: (energy eV, weight)} -> {key: population}, normalised"""
    e0 = min(e for e, g in levels.values())
    p = {k: g * math.exp(-(e - e0) / (L.KB_EV * T)) for k, (e, g) in levels.items()}
    s = sum(p.values())
    return {k: v / s for k, v in p.items()}


def build(gid, Tg, inp):
    files, mkey = SETS[gid]
    blocks = []
    for f in files:
        blocks += L.read_file(os.path.join(inp, f))
    # LoKI-MC's Collision::add keeps the first of two identical collisions
    seen, colls = set(), []
    for b in blocks:
        key = (b['type'], b['target'], tuple(sorted(b['products'])), b['thr'], b['rev'], b['sub'])
        if key in seen:
            continue
        seen.add(key)
        colls.append(b)
    M = L.db_value(inp, 'masses.txt', mkey) if mkey else 4.000000 * 1.660539040e-27
    names = set()
    for b in colls:
        names.add(b['target'])
        for c, n in b['products']:
            names.add(n)
    ground = [b['target'] for b in colls if b['type'] in ('Elastic', 'Effective')][0]

    def pops(T):
        P = {n: 0.0 for n in names}
        G = {n: None for n in names}
        P[ground] = 1.0
        G[ground] = 1.0
        g, ele, _, _ = state_parts(ground)
        # the vibrational children of the ground electronic state, and the rotational children of its v = 0 (a
        # vibrational label on an excited electronic state, N2(A3Su+,v=10+), is not one)
        vib = {n: state_parts(n)[2] for n in names if state_parts(n)[1] == ele and state_parts(n)[2] is not None
               and state_parts(n)[3] is None}
        rot = {n: state_parts(n)[3] for n in names if state_parts(n)[1] == ele and state_parts(n)[3] is not None}
        if gid == 'CO2':
            P.update({'CO2(X)': 1.0, 'CO2(X,v=000)': 0.9254, 'CO2(X,v=010)': 0.0746})
            for n, w in [('CO2(X,v=000)', 1), ('CO2(X,v=010)', 2), ('CO2(X,v=020)', 3), ('CO2(X,v=100)', 1),
                         ('CO2(X,v=001)', 1), ('CO2(X,v=030+110)', 6), ('CO2(X,v=040+120+011)', 10), ('CO2(X,v=200)', 1),
                         ('CO2(X,v=050+210+130+021+101)', 16), ('CO2(X,v=300)', 1), ('CO2(X,v=060+220+140)', 15)]:
                G[n] = w
            return P, G
        if vib and gid in ('N2', 'O2', 'CO', 'H2O'):
            lev = {}
            for n, v in vib.items():
                if gid == 'N2':
                    e = {'0': 0.0}.get(v)
                    if e is None:
                        e = N2_VIB[v]
                    lev[n] = (e, 1.0)
                elif gid == 'O2':
                    w = L.db_value(inp, 'harmonicFrequencies.txt', 'O2')
                    wx = L.db_value(inp, 'anharmonicFrequencies.txt', 'O2')
                    vv = int(v) + 0.5
                    lev[n] = (HBAR_EV * (w * vv - wx * vv * vv), 3.0)
                elif gid == 'CO':
                    w = L.db_value(inp, 'harmonicFrequencies.txt', 'CO')
                    lev[n] = (HBAR_EV * w * (int(v) + 0.5), 1.0)
                elif gid == 'H2O':
                    lev[n] = (H2O_VIB[v], 1.0)
            for n, p in boltz(lev, T).items():
                P[n] = p
                G[n] = lev[n][1]
        if rot:
            B = L.db_value(inp, 'rotationalConstants.txt', gid)
            lev = {}
            for n, J in rot.items():
                J = int(J)
                g = 3.0 * (1 + 0.5 * (1 + (-1) ** J)) * (2 * J + 1) if gid == 'N2' else 2.0 * J + 1
                lev[n] = (B * J * (J + 1), g)
            v0 = '%s(X,v=0)' % gid
            pv0 = P.get(v0, 1.0)
            for n, p in boltz(lev, T).items():
                P[n] = pv0 * p
                G[n] = lev[n][1]
        return P, G

    P, Gw = pops(Tg)
    P300, G300 = pops(300.0)
    eff = [b for b in colls if b['type'] == 'Effective']
    out = []
    if eff:
        eff = eff[0]
        x = eff['e']
        s = eff['s'].copy()
        for b in colls:
            if b['type'] in ('Effective', 'Elastic'):
                continue
            s -= P300.get(b['target'], 0.0) * L.interp(b, x)
            if b['rev'] and P300.get(b['products'][0][1], 0.0) > 0:
                gr = G300[b['target']] / G300[b['products'][0][1]]
                s -= P300[b['products'][0][1]] * L.superelastic(b, x, gr)
        neg = int(np.sum(s < 0))
        s = np.where(s < 0, 0.0, s)
        el = dict(eff)
        el.update(type='Elastic', e=x, s=s, thr=0.0, derived=True, clipped=neg)
    else:
        el = [b for b in colls if b['type'] == 'Elastic' and P.get(b['target'], 0) > 0][0]
    for b in colls:
        if b['type'] in ('Effective', 'Elastic'):
            continue
        c = dict(b)
        c['w'] = P.get(b['target'], 0.0)
        if b['rev']:
            prod = b['products'][0][1]
            c['wsup'] = P.get(prod, 0.0)
            if c['wsup'] > 0:
                c['gratio'] = Gw[b['target']] / Gw[prod]
        out.append(c)
    return dict(id=gid, M=M, el=el, colls=out, sigma_el=lambda e, el=el: L.interp(el, e) * P.get(el['target'], 1.0),
                pops=P, ground=ground, n_blocks=len(blocks), n_colls=len(colls))


N2_VIB = {}
H2O_VIB = {}


def load_tables(inp):
    for ln in open(os.path.join(inp, 'Nitrogen', 'N2energyVibLevels.txt')):
        t = ln.split()
        if len(t) == 2 and t[0].startswith('N2(X,v='):
            N2_VIB[t[0][7:-1]] = float(t[1])
    for ln in open(os.path.join(inp, 'H2O', 'H2O_vibEnergy_Itikawa.txt')):
        t = ln.split()
        if len(t) == 2 and t[0].startswith('H2O(X,v='):
            H2O_VIB[t[0][8:-1]] = float(t[1])
