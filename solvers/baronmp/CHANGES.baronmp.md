Summary of recent updates to BARONMP for AMPL
=============================================


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).


## 20250811
- Fixed floating-point output precision for model
  submission to Baron.
- Changes in MP
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


## 20250806
- Updated Baron to version 2025.8.5, which includes
  bugfixes.


## 20250801
- Changes in MP
  - Tolerances set by options pre:feastol,
    pre:feastolrel both need to be violated
    to produce a warning on contradicting
    variable/constraint bounds. Previously
    the preprocessor failed on any violation,
    without letting the solver try.
  - Options cvt:compl, cvt:compl:eps control
    complementarity reformulations.


## 20250617
- Changes in MP
  - Multi-objective emulator: added support for
    objective-specific options via objective suffixes
    beginning with *option_*
  - Option *cvt:unnest*: bits 2 and 4 switch on
    inlining of linear and quadratic subexpressions
    produced during reformulations (by default on).
  - Options *cvt:pre:ctx2ineq*, *cvt:pre:ctx2count*
    to control context propagation into conditional
    comparisons #267.



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
