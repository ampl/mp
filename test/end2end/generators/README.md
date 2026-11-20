# Generators of test models

`logical_models.py` generates the random logical models in
`../cases/categorized/fast_random`: small integer programs whose expected results come from
enumerating all integer points of their bounded domains. A round is determined by its
generator and number, so existing rounds are reproduced exactly, and new rounds give new
models.

| Generator | Models per round | Seed | Content |
| --- | --- | --- | --- |
| `glog` | 200 | 3000 + round | implications, equivalences, `and`/`or`/`not`, `alldiff`, `numberof`, `atleast`/`atmost`, `forall`/`exists`, `if`-`then`-`else`, `max`/`min`, `abs` |
| `glog2` | 500 | 5000 + round | the `glog` constraints in ten categories (model number modulo 10): default, constants at the domain edges, repeated variables, fixed variables, wide negative domains, many constraints, feasibility only, large coefficients, binary, duplicated constraints |
| `gctx` | 200 | 4000 + round | the same atoms in positive, negative, mixed and nested contexts, and in the objective |

Generate a round (from this folder; it writes `<generator>_r<round>/` with `m_<k>.mod`
and `modellist.json`):

```
python logical_models.py glog 6
python logical_models.py gctx 4 --count 50 --out /tmp/cases
```

Rounds 1-5 of `glog` and `glog2` and rounds 1-3 of `gctx` are in `fast_random`. Run a
round with the end-to-end runner, for example from `test/end2end`:

```
python run.py highs --dir cases/categorized/fast_random/glog_r6 --nonrecursive
```

Layout: the models are written with `amplfmt.py` (at most 80 characters per line, long
statements broken before the operator of lowest precedence). `--raw` writes one statement
per line, as the rounds were first committed. `nlcheck.py` checks that a layout change
leaves the NL files unchanged; it needs AMPL (`AMPL` environment variable or `--ampl`):

```
python nlcheck.py ../cases/categorized/fast_random/glog_r6
python nlcheck.py --careful <folder of third-party models>
```

Enumeration takes up to a few seconds per model (up to 7^6 points); a round of 500
models takes a few minutes.
