# Randomly generated logical models

Small integer models with logical constraints and expressions, generated at random to test
the reformulation of logical constraints in MP. They are not run by default: point the
end-to-end runner at this folder, or at one of its subfolders, with `--dir`. For example,
from `test/end2end`:

```
python run.py highs --dir cases/categorized/fast_random/glog_r1 --nonrecursive
```

Each subfolder is one generated round, with its own `modellist.json`:

| Folders | Models | What they exercise |
| --- | --- | --- |
| `glog_r1` ... `glog_r5` | 200 each | implications, equivalences, `and`/`or`/`not`, `alldiff`, `numberof`, `atleast`/`atmost`, comparisons in the objective |
| `glog2_r1` ... `glog2_r5` | 500 each | special categories: boundary constants, repeated variables, fixed variables, wide negative domains, many constraints, feasibility only, large coefficients, nested expressions |
| `gctx_r1` ... `gctx_r3` | 200 each | the same atoms (comparisons, `abs`, `max`, `if`-`then`-`else`) in positive, negative and mixed contexts, nested, and in the objective |

The expected results come from enumerating all points of the bounded integer domain. A
case expects either the optimal objective value or `solve_result` `infeasible`. A failure
means that the driver, the conversion or the solver returned something else; a solver's
own tolerances can matter for large coefficients.

`../../../generators/logical_models.py` generates the rounds (and new ones); see the
README there. The models were first committed as generated, one statement per line, and
then laid out to be read (at most 80 characters per line, blank lines between
declarations, objectives and constraints). The layout change leaves the NL file of every
model unchanged.
