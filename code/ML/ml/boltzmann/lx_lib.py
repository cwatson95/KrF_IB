#!/usr/bin/env python3
"""lx_lib.py -- derive.py's own reading of LoKI-MC's cross sections, written apart from extract.py and from the page.

It reads a file as LoKI-MC's Parse::LXCatFiles does (the block opens at a 'PROCESS:' line; the threshold is the
'E = <x> eV' of the next line; the description LoKI-MC parses is the line after; the table sits between the next two
'-----' lines) but by a regular expression over the whole text rather than a line loop, and it parses the bracketed
description as Parse::addLXCatEntry does (spaces removed, the type the last comma field, '->' or '<->', states split at
'+' outside parentheses, a leading integer a stoichiometric coefficient).  The interpolation is LoKI-MC's
Collision::interpolatedCrossSection: linear, zero outside the table, and zero at or below the threshold of any
collision that is not Elastic or Effective.  An Effective cross section becomes an Elastic one as
EedfGas::elasticFromEffectiveCrossSection does it, on the effective table's own energies, clipped at zero.

Constants are LoKI-MC's (Headers/Constant.h, CODATA 2014), so that its rate coefficients are reproduced to rounding.
"""
import math, os, re
import numpy as np

QE = 1.6021766208e-19      # LoKI-MC Constant::electronCharge
ME = 9.10938356e-31        # Constant::electronMass
KB = 1.38064852e-23        # Constant::boltzmann
AMU = 1.660539040e-27      # Constant::unifiedAtomicMass
KB_EV = KB / QE
GAMMA = math.sqrt(2.0 * QE / ME)   # sqrt(2e/m): v = GAMMA * sqrt(eps[eV])
TYPES = ('Elastic', 'Effective', 'Attachment', 'Ionization', 'Excitation', 'Vibrational', 'Rotational')

BLOCK_RE = re.compile(r'^(.*PROCESS:.*)\n(.*)\n(.*)\n(?:.*\n)*?.*-----.*\n((?:(?!.*-----).*\n)*)', re.M)


def split_states(s):
    """split 'A(x,y)+B+2C(z)' at '+' outside parentheses; returns [(coeff, name), ...]"""
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == '+' and depth == 0 and cur:
            out.append(cur)
            cur = ''
        else:
            cur += ch
    if cur:
        out.append(cur)
    res = []
    for t in out:
        m = re.match(r'^(\d+)(?=[A-Za-z])', t)
        c = int(m.group(1)) if m else 1
        res.append((c, t[m.end():] if m else t))
    return res


def parse_desc(desc):
    d = desc.replace(' ', '')
    body = d[d.index('[') + 1:d.rindex(']')]
    fields = body.split(',')
    sub = None
    if fields[-1].lower() in ('momentum-transfer', 'integral'):
        sub = fields.pop().lower()
    typ = fields.pop()
    assert typ in TYPES, 'unknown collision type in ' + desc
    reac = ','.join(fields)
    rev = '<->' in reac
    lhs, rhs = reac.split('<->' if rev else '->')
    L = [x for x in split_states(lhs) if x[1] not in ('e', 'E')]
    R = split_states(rhs)
    ne = sum(c for c, n in R if n in ('e', 'E'))
    R = [x for x in R if x[1] not in ('e', 'E')]
    assert len(L) == 1, 'one target expected: ' + desc
    if sub is None:
        sub = 'momentum-transfer' if typ in ('Elastic', 'Effective') else 'integral'
    return dict(type=typ, sub=sub, target=L[0][1], products=R, rev=rev, eout=ne, reaction=reac)


def read_file(path):
    text = open(path, 'rb').read().decode('utf-8', errors='replace')
    if not text.endswith('\n'):
        text += '\n'
    blocks = []
    for m in BLOCK_RE.finditer(text):
        proc, param, desc, tab = m.group(1), m.group(2), m.group(3), m.group(4)
        t = re.search(r'E =(.*?)eV', param)
        thr = float(t.group(1)) if t else 0.0
        rows = [ln.split() for ln in tab.splitlines() if ln.strip()]
        e = np.array([float(r[0]) for r in rows])
        s = np.array([float(r[1]) for r in rows])
        b = dict(proc=proc.strip(), param=param.strip(), desc=desc.strip(), thr=thr, e=e, s=s)
        b.update(parse_desc(desc))
        blocks.append(b)
    return blocks


def interp(b, x):
    """LoKI-MC Collision::interpolatedCrossSection, the 'integral' component (= the table; isotropic)"""
    x = np.asarray(x, dtype=float)
    e, s = b['e'], b['s']
    y = np.interp(x, e, s)
    y = np.where((x >= e[0]) & (x <= e[-1]), y, 0.0)
    if b['type'] not in ('Elastic', 'Effective'):
        y = np.where(x > b['thr'], y, 0.0)
    return y


def superelastic(b, x, g_ratio):
    """Collision::superElasticCrossSection: (g_target/g_product) (1 + u/x) sigma(x + u), zero at x = 0"""
    x = np.asarray(x, dtype=float)
    u = b['thr']
    y = np.zeros_like(x)
    nz = x > 0
    y[nz] = g_ratio * (1.0 + u / x[nz]) * interp(b, x[nz] + u)
    return y


def eval_expr(s):
    """the Databases/ expressions ('39.948*1.660539040e-27', '2*pi*...'): numbers, * / + - ( ) and pi only"""
    assert re.fullmatch(r'[0-9eE.+\-*/() pi]+', s), 'unexpected expression ' + s
    return eval(s, {'__builtins__': {}}, {'pi': math.pi})


def db_value(inp, fname, key):
    for ln in open(os.path.join(inp, 'Databases', fname), encoding='utf-8'):
        t = ln.split('%')[0].split()
        if len(t) >= 2 and t[0] == key:
            return eval_expr(''.join(t[1:]))
    raise KeyError(key)


# ---------------- the rate coefficient as LoKI-MC evaluates it from an EEDF on its grid ----------------
def loki_rate(b, cells, f, sigma_nodes=None):
    """Collision::evaluateRateCoeff: k = sqrt(2e/m) sum_c 0.5(s(u_c-) + s(u_c+)) u_c f_c du over the cells from
    floor(threshold/du) up; cells are the cell centres of eedf.txt, the nodes 0, du, 2du, ..."""
    du = 2.0 * cells[0]
    nodes = np.arange(len(cells) + 1) * du
    s = interp(b, nodes) if sigma_nodes is None else sigma_nodes
    lmin = int(math.floor(b['thr'] / du))
    if len(cells) - lmin <= 0:
        return 0.0
    cs = 0.5 * (s[lmin:-1] + s[lmin + 1:])
    return GAMMA * float(np.sum(cs * cells[lmin:] * f[lmin:])) * du


# ---------------- LoKI-MC output files ----------------
def read_eedf(path):
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.lstrip().startswith('Energy')]
    a = np.array([[float(x) for x in r[:2]] for r in rows])
    return a[:, 0], a[:, 1]


def read_rates(path):
    out = []
    for ln in open(path):
        t = ln.split()
        if not t or not t[0].isdigit():
            continue
        ine = float(t[1])
        sup = None if t[2] == '(N/A)' else float(t[2])
        out.append(dict(id=int(t[0]), ine=ine, sup=sup, desc=' '.join(t[3:])))
    return out


def read_swarm(path):
    txt = open(path).read()
    def grab(label, which=0):
        ms = re.findall(re.escape(label) + r'\s*=\s*([-+0-9.eE]+)', txt)
        return float(ms[which]) if len(ms) > which else None
    vz = re.findall(r'\| v_z\' \|\s*\|\s*([-+0-9.eE]+)', txt)
    return dict(EN=grab('Reduced electric field'), mean=grab('Mean energy'), Te=grab('Electron temperature'),
                mu_flux=grab('Reduced mobility coefficient', 0), mu_bulk=grab('Reduced mobility coefficient', 1),
                DT_flux=grab('Reduced transverse diffusion coefficient', 0),
                char_flux=grab('Characteristic energy', 0), alpha_flux=grab('Reduced Townsend coefficient', 0),
                alpha_bulk=grab('Reduced Townsend coefficient', 1), eta_flux=grab('Reduced attachment coefficient', 0),
                kion=grab('Ionization coefficient'), katt=grab('Attachment coefficient'),
                mu_eedf=grab('DC non-magnetized reduced mobility'), DN_eedf=grab('Reduced diffusion coefficient'),
                char_eedf=grab('Characteristic energy', 2),
                vz_flux=abs(float(vz[0])) if vz else None, vz_bulk=abs(float(vz[1])) if len(vz) > 1 else None)
