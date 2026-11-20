#!/usr/bin/env python3
"""Generate the random logical models of cases/categorized/fast_random.

Each model is a small integer program over x{1..n} with a few domain values; its
expected result comes from enumerating all integer points (optimal objective value,
or infeasible). Three generators, each seeded by the round number, so a round is
reproduced exactly:

  glog   200 models per round, seed 3000 + R: implications, equivalences, and/or/not,
         alldiff, numberof, atleast/atmost, forall/exists, if-then-else, max/min, abs
  glog2  500 models per round, seed 5000 + R: the glog constraints in special
         categories (model number modulo 10): 0 default, 1 constants at the domain
         edges, 2 the same variable repeated, 3 fixed variables, 4 wide negative
         domains, 5 many constraints, 6 feasibility only, 7 large coefficients,
         8 binary, 9 duplicated constraints
  gctx   200 models per round, seed 4000 + R: the same atoms (comparisons, abs,
         max/min, if-then-else) in positive, negative, mixed and nested contexts,
         and in the objective

Usage (from test/end2end/generators):

  python logical_models.py glog 6               -> ../cases/categorized/fast_random/glog_r6
  python logical_models.py gctx 4 --count 50 --out /tmp/cases
  python logical_models.py glog2 1 --raw        one statement per line (as first committed)

The output folder gets m_<k>.mod and a modellist.json for run.py. By default the
models are laid out with amplfmt.py, which changes the layout only.
"""
import argparse
import itertools
import json
import os
import random
import re

from amplfmt import format_model

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(HERE, "..", "cases", "categorized", "fast_random")


# ---- glog and glog2: one random constraint of a template --------------------------

# glog2 sets these; "v0": the extra draw added for glog2's repeated variables (the glog
# rounds were generated before it, so they are reproduced without it)
MODE = {"same_index": False, "edge_consts": False, "v0": False}   # set by glog2 for special categories

def templ(rnd, n, lo, hi):
    """one random constraint: (ampl text, predicate over a list x of values)"""
    _v0 = rnd.randrange(n) if MODE["v0"] else 0
    V = (lambda: _v0 if MODE["same_index"] and rnd.random() < 0.6 else rnd.randrange(n))
    C = (lambda: rnd.choice([lo, hi, lo - 1, hi + 1]) if MODE["edge_consts"] else rnd.randint(lo, hi))
    op = lambda: rnd.choice(["<=", ">=", "==", "!=", "<", ">"])
    PY = {"<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b, "==": lambda a, b: a == b,
          "!=": lambda a, b: a != b, "<": lambda a, b: a < b, ">": lambda a, b: a > b}
    xs = lambda i: f"x[{i+1}]"
    k = rnd.randrange(16)
    if k == 0:
        a, b, o1, o2, c1, c2 = V(), V(), op(), op(), C(), C()
        return f"({xs(a)} {o1} {c1}) ==> ({xs(b)} {o2} {c2})", lambda x: (not PY[o1](x[a], c1)) or PY[o2](x[b], c2)
    if k == 1:
        a, b, c, o1, o2, c1, c2 = V(), V(), V(), op(), op(), C(), C()
        return f"({xs(a)} {o1} {c1}) <==> ({xs(b)} + {xs(c)} {o2} {c2})", lambda x: PY[o1](x[a], c1) == PY[o2](x[b] + x[c], c2)
    if k == 2:
        a, b, c, o1, o2, o3, c1, c2, c3 = V(), V(), V(), op(), op(), op(), C(), C(), C()
        return f"({xs(a)} {o1} {c1}) or ({xs(b)} {o2} {c2}) or ({xs(c)} {o3} {c3})", lambda x: PY[o1](x[a], c1) or PY[o2](x[b], c2) or PY[o3](x[c], c3)
    if k == 3:
        a, b, o1, o2, c1, c2 = V(), V(), op(), op(), C(), C()
        return f"not (({xs(a)} {o1} {c1}) and ({xs(b)} {o2} {c2}))", lambda x: not (PY[o1](x[a], c1) and PY[o2](x[b], c2))
    if k == 4:
        a, b, c, o1, o2, o3, c1, c2, c3 = V(), V(), V(), op(), op(), op(), C(), C(), C()
        return (f"({xs(a)} {o1} {c1}) ==> ({xs(b)} {o2} {c2}) else ({xs(c)} {o3} {c3})",
                lambda x: PY[o2](x[b], c2) if PY[o1](x[a], c1) else PY[o3](x[c], c3))
    if k == 5:
        kk, cnt, o = C(), rnd.randint(0, 3), op()
        return f"numberof {kk} in ({', '.join(xs(i) for i in range(n))}) {o} {cnt}", lambda x: PY[o](sum(1 for v in x if v == kk), cnt)
    if k == 6:
        s = rnd.sample(range(n), rnd.randint(2, min(4, n)))
        return f"alldiff({', '.join(xs(i) for i in s)})", lambda x: len({x[i] for i in s}) == len(s)
    if k == 7:
        t, kk = C(), rnd.randint(1, 3)
        if rnd.random() < 0.5:
            return f"atleast {kk} {{i in 1..{n}}} (x[i] >= {t})", lambda x: sum(1 for v in x if v >= t) >= kk
        return f"atmost {kk} {{i in 1..{n}}} (x[i] >= {t})", lambda x: sum(1 for v in x if v >= t) <= kk
    if k == 8:
        return f"forall {{i in 1..{n-1}}} x[i] != x[i+1]", lambda x: all(x[i] != x[i + 1] for i in range(n - 1))
    if k == 9:
        t = C()
        return f"exists {{i in 1..{n}}} x[i] == {t}", lambda x: any(v == t for v in x)
    if k == 10:
        a, b, c, o1, o2, c1, c2 = V(), V(), V(), op(), op(), C(), C()
        return f"(if {xs(a)} {o1} {c1} then {xs(b)} else {xs(c)}) {o2} {c2}", lambda x: PY[o2](x[b] if PY[o1](x[a], c1) else x[c], c2)
    if k == 11:
        a, b, c, d, e, o = V(), V(), V(), V(), C(), op()
        return f"max({xs(a)}, {xs(b)}) - min({xs(c)}, {xs(d)}) {o} {e}", lambda x: PY[o](max(x[a], x[b]) - min(x[c], x[d]), e)
    if k == 12:
        a, b, d, o = V(), V(), rnd.randint(0, hi - lo), op()
        return f"abs({xs(a)} - {xs(b)}) {o} {d}", lambda x: PY[o](abs(x[a] - x[b]), d)
    if k == 13:
        a, b, c, o1, o2, c1, c2 = V(), V(), V(), op(), op(), C(), C()
        return (f"({xs(a)} {o1} {c1}) and ({xs(b)} {o2} {c2}) ==> ({xs(c)} >= {c1})",
                lambda x: (not (PY[o1](x[a], c1) and PY[o2](x[b], c2))) or x[c] >= c1)
    if k == 14:
        a, b, o1, c1 = V(), V(), op(), C()
        return f"((if {xs(a)} {o1} {c1} then 1 else 0) + (if {xs(b)} >= {c1} then 1 else 0)) >= 1", lambda x: (PY[o1](x[a], c1) + (x[b] >= c1)) >= 1
    a, b, c, o1, o2, c1, c2 = V(), V(), V(), op(), op(), C(), C()
    return (f"(({xs(a)} {o1} {c1}) ==> ({xs(b)} {o2} {c2})) and (({xs(b)} {o2} {c2}) ==> ({xs(c)} >= {c1}))",
            lambda x: ((not PY[o1](x[a], c1)) or PY[o2](x[b], c2)) and ((not PY[o2](x[b], c2)) or x[c] >= c1))


# ---- gctx: atoms reused in different contexts -------------------------------------

OPS = {"<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b, "==": lambda a, b: a == b, "!=": lambda a, b: a != b,
       "<": lambda a, b: a < b, ">": lambda a, b: a > b}

def atom(rnd, n, lo, hi):
    """a Boolean atom: (ampl text, predicate)"""
    V = lambda: rnd.randrange(n); C = lambda: rnd.randint(lo, hi); op = lambda: rnd.choice(list(OPS))
    k = rnd.randrange(7); a, b, c, o = V(), V(), C(), op()
    if k < 2: return f"x[{a+1}] {o} {c}", lambda x: OPS[o](x[a], c)
    if k == 2: return f"x[{a+1}] + x[{b+1}] {o} {c}", lambda x: OPS[o](x[a] + x[b], c)
    if k == 3: return f"abs(x[{a+1}] - x[{b+1}]) {o} {rnd.randint(0, hi - lo)}", None
    if k == 4: return f"max(x[{a+1}], x[{b+1}]) {o} {c}", lambda x: OPS[o](max(x[a], x[b]), c)
    if k == 5: return f"min(x[{a+1}], x[{b+1}]) {o} {c}", lambda x: OPS[o](min(x[a], x[b]), c)
    return f"x[{a+1}] {o} x[{b+1}]", lambda x: OPS[o](x[a], x[b])
def mk_atoms(rnd, n, lo, hi, k):
    out = []
    for _ in range(k):
        while True:
            t, p = atom(rnd, n, lo, hi)
            if p is None:                     # abs(x[a]-x[b]) o d  (parse back)
                m = re.match(r"abs\(x\[(\d+)\] - x\[(\d+)\]\) (\S+) (\d+)", t); a, b, o, d = int(m[1]) - 1, int(m[2]) - 1, m[3], int(m[4])
                p = (lambda a, b, o, d: lambda x: OPS[o](abs(x[a] - x[b]), d))(a, b, o, d)
            out.append((t, p)); break
    return out

def constraints(rnd, atoms, n, lo, hi):
    """contexts: positive, negative, mixed, nested. Atoms are reused verbatim to provoke shared flat expressions."""
    A = lambda: rnd.choice(atoms); V = lambda: rnd.randrange(n); C = lambda: rnd.randint(lo, hi)
    cons = []
    for _ in range(rnd.randint(2, 4)):
        k = rnd.randrange(12); a1, a2, a3 = A(), A(), A()
        P = lambda a: a[1]
        if k == 0: cons.append((f"({a1[0]}) or ({a2[0]})", lambda x, a1=a1, a2=a2: a1[1](x) or a2[1](x)))
        elif k == 1: cons.append((f"not ({a1[0]})", lambda x, a1=a1: not a1[1](x)))
        elif k == 2: cons.append((f"({a1[0]}) ==> ({a2[0]})", lambda x, a1=a1, a2=a2: (not a1[1](x)) or a2[1](x)))
        elif k == 3: cons.append((f"({a1[0]}) <==> ({a2[0]})", lambda x, a1=a1, a2=a2: a1[1](x) == a2[1](x)))
        elif k == 4:
            i, j, c, o = V(), V(), C(), rnd.choice(["<=", ">=", "<", ">", "=="])
            cons.append((f"(if {a1[0]} then x[{i+1}] else x[{j+1}]) {o} {c}", lambda x, a1=a1, i=i, j=j, c=c, o=o: OPS[o](x[i] if a1[1](x) else x[j], c)))
        elif k == 5:
            c = rnd.randint(1, 3)
            cons.append((f"(if {a1[0]} then 1 else 0) + (if {a2[0]} then 1 else 0) + (if {a3[0]} then 1 else 0) >= {c}",
                         lambda x, a1=a1, a2=a2, a3=a3, c=c: (a1[1](x) + a2[1](x) + a3[1](x)) >= c))
        elif k == 6: cons.append((f"not (({a1[0]}) and ({a2[0]}))", lambda x, a1=a1, a2=a2: not (a1[1](x) and a2[1](x))))
        elif k == 7: cons.append((f"({a1[0]}) ==> (({a2[0]}) or ({a3[0]}))", lambda x, a1=a1, a2=a2, a3=a3: (not a1[1](x)) or a2[1](x) or a3[1](x)))
        elif k == 8: cons.append((f"({a1[0]}) ==> ({a2[0]}) else ({a3[0]})", lambda x, a1=a1, a2=a2, a3=a3: a2[1](x) if a1[1](x) else a3[1](x)))
        elif k == 9: cons.append((f"({a1[0]}) and (({a2[0]}) or not ({a3[0]}))", lambda x, a1=a1, a2=a2, a3=a3: a1[1](x) and (a2[1](x) or not a3[1](x))))
        elif k == 10:
            i, j, o, c = V(), V(), rnd.choice(["<=", ">=", "<", ">"]), C()
            cons.append((f"(if {a1[0]} then x[{i+1}] else 0) + (if {a2[0]} then x[{j+1}] else 0) {o} {c}",
                         lambda x, a1=a1, a2=a2, i=i, j=j, o=o, c=c: OPS[o]((x[i] if a1[1](x) else 0) + (x[j] if a2[1](x) else 0), c)))
        else: cons.append((f"({a1[0]}) <==> ({a2[0]})", lambda x, a1=a1, a2=a2: a1[1](x) == a2[1](x)))
    return cons


# ---- the generators: (AMPL text, expected result) per model ------------------------


def best_value(n, lo, hi, cons, objective, sense):
    """Optimal objective value over the integer points of [lo, hi]^n, or None."""
    best = None
    for x in itertools.product(range(lo, hi + 1), repeat=n):
        if all(p(x) for _, p in cons):
            v = objective(x)
            if best is None or (v < best if sense == "minimize" else v > best):
                best = v
    return best


def gen_glog(R, N=200):
    rnd = random.Random(3000 + R)
    for _ in range(N):
        n = rnd.randint(3, 6); lo = rnd.choice([0, 0, -2]); hi = lo + rnd.randint(2, 4)
        cons = [templ(rnd, n, lo, hi) for _ in range(rnd.randint(1, 4))]
        coef = [rnd.randint(-5, 6) for _ in range(n)]; sense = rnd.choice(["minimize", "maximize"])
        t = f"var x{{1..{n}}} integer >= {lo} <= {hi};\n"
        t += f"{sense} o: " + " + ".join(f"({c})*x[{i+1}]" for i, c in enumerate(coef)) + ";\n"
        for j, (s, _) in enumerate(cons):
            t += f"s.t. c{j+1}: {s};\n"
        yield t, best_value(n, lo, hi, cons, lambda x: sum(c * xi for c, xi in zip(coef, x)), sense)


def gen_glog2(R, N=500):
    MODE["v0"] = True
    try:
        yield from _gen_glog2(R, N)
    finally:
        MODE["v0"] = False


def _gen_glog2(R, N):
    rnd = random.Random(5000 + R)
    for k in range(N):
        cat = k % 10
        n = rnd.randint(3, 6); lo = rnd.choice([0, 0, -2]); hi = lo + rnd.randint(2, 4)
        if cat == 4:
            lo = -rnd.randint(2, 3); hi = lo + rnd.randint(4, 5); n = min(n, 5)
        if cat == 8:
            lo, hi, n = 0, 1, rnd.randint(4, 7)
        nc = rnd.randint(5, 7) if cat == 5 else rnd.randint(1, 4)
        MODE["same_index"] = (cat == 2); MODE["edge_consts"] = (cat == 1)
        cons = [templ(rnd, n, lo, hi) for _ in range(nc)]
        MODE["same_index"] = MODE["edge_consts"] = False
        if cat == 9 and cons:
            cons = cons + [cons[0]] + ([cons[-1]] if len(cons) > 1 else [])
        if cat == 3:
            for i in rnd.sample(range(n), rnd.randint(1, 2)):
                c = rnd.randint(lo, hi)
                cons.append((f"x[{i+1}] == {c}", (lambda i, c: lambda x: x[i] == c)(i, c)))
        coef = [rnd.randint(-5, 6) for _ in range(n)]
        if cat == 7:
            coef = [rnd.choice([-1000, 1000, 17, -333]) * rnd.randint(1, 3) for _ in range(n)]
        if cat == 6:
            coef = [0] * n
        sense = rnd.choice(["minimize", "maximize"])
        t = (f"var x{{1..{n}}} integer >= {lo} <= {hi};\n{sense} o: "
             + " + ".join(f"({c})*x[{i+1}]" for i, c in enumerate(coef)) + ";\n")
        for j, (s, _) in enumerate(cons):
            t += f"s.t. c{j+1}: {s};\n"
        yield t, best_value(n, lo, hi, cons, lambda x: sum(c * xi for c, xi in zip(coef, x)), sense)


def gen_gctx(R, N=200):
    rnd = random.Random(4000 + R)
    for _ in range(N):
        n = rnd.randint(3, 5); lo = rnd.choice([0, 0, -1]); hi = lo + rnd.randint(2, 3)
        atoms = mk_atoms(rnd, n, lo, hi, rnd.randint(2, 4))
        cons = constraints(rnd, atoms, n, lo, hi)
        coef = [rnd.randint(-5, 6) for _ in range(n)]; sense = rnd.choice(["minimize", "maximize"])
        obj = " + ".join(f"({c})*x[{i+1}]" for i, c in enumerate(coef))
        extra = None
        if rnd.random() < 0.3:                       # an atom in the objective too
            a = rnd.choice(atoms); w = rnd.randint(-6, 6)
            obj += f" + ({w})*(if {a[0]} then 1 else 0)"; extra = (a, w)
        t = (f"var x{{1..{n}}} integer >= {lo} <= {hi};\n{sense} o: {obj};\n"
             + "".join(f"s.t. c{j+1}: {s};\n" for j, (s, _) in enumerate(cons)))

        def objective(x, coef=coef, extra=extra):
            v = sum(c * xi for c, xi in zip(coef, x))
            return v + (extra[1] * (1 if extra[0][1](x) else 0) if extra else 0)
        yield t, best_value(n, lo, hi, cons, objective, sense)


GENERATORS = {"glog": (gen_glog, 200), "glog2": (gen_glog2, 500), "gctx": (gen_gctx, 200)}


def write_round(kind, R, count, out, raw):
    gen, default = GENERATORS[kind]
    name = f"{kind}_r{R}"
    folder = os.path.join(out, name)
    os.makedirs(folder, exist_ok=True)
    cases = []
    for k, (text, best) in enumerate(gen(R, count or default)):
        with open(os.path.join(folder, f"m_{k}.mod"), "w") as f:
            f.write(text if raw else format_model(text))
        case = {"name": f"{name}_m_{k}", "files": [f"m_{k}.mod"], "tags": ["logical", "integer"]}
        if best is None:
            case["values"] = {"solve_result": "infeasible"}
        else:
            case["objective"] = best
        cases.append(case)
    with open(os.path.join(folder, "modellist.json"), "w") as f:
        json.dump(cases, f, indent=1)
    print(f"{folder}: {len(cases)} models")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("kind", choices=sorted(GENERATORS))
    ap.add_argument("round", type=int)
    ap.add_argument("--count", type=int, help="models (default 200, glog2 500)")
    ap.add_argument("--out", default=DEFAULT_OUT, help="parent folder of <kind>_r<round>")
    ap.add_argument("--raw", action="store_true", help="no layout: one statement per line")
    a = ap.parse_args()
    write_round(a.kind, a.round, a.count, a.out, a.raw)


if __name__ == "__main__":
    main()
