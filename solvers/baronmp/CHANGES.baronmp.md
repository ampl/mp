Summary of recent updates to BARONMP for AMPL
=============================================


## unreleased
- Changes in MP
  - Option cvt:qp2pass (default even faster parsing
    of quadratics)
  - Fix a bug in parsing of quadratic expressions,
    which could wrongly parse products of unequal
    linear expressions, such as (x-3)*(x-z-5).


## 20250329
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.


## 20250308
- Changes in MP.


## 20241227
- Updated to baron 24.12.21. In addition to bug fixes, this version comes with improvements in memory management, convexity identification, branching, and new relaxation and reduction strategies for certain quadratic and concave programs.


## 20241209
- Removed spurious output from execution of local solves.
- Fixed bug that made `outlev` ineffective.


## 20241206
- Fixed issue with variable and constraint names containing invalid
  characters
- Added solve results information


## 20241119
- Beta release of MP driver for Baron, use with option `option solver baronmp;`,
  to use the original ASL driver use `option solver baron;`
- Breaking changes with ASL driver:
  - `version` keyword removed
  - Changed behaviour of the `iisfind` option to be consistent with MP features;
    see `alg:iisfind`, `alg:iismethod`, `alg:iisint`
  - Absolute values are reformulated as MIP instead of abs(x)=(x^2)^0.5
  - Added option value `local` to `scratch` that automatically places the scratch 
    files in a directory with the name derived from the model file (e.g. for 
    `/test/model.nl` the directory will be `/test/model.nl_baron/`)
