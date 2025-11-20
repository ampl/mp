Summary of recent updates to cuOpt for AMPL
===========================================


## 20251121
- New option *obj:multi:options* to control
  whether multi-objective option suffixes
  are used.


## 20251119
- Updated to cuOpt 25.10.1
  - Barrier algorithm (see option *alg:method*)
  
  
## 20251015
- Changes in MP


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options *cvt:pre:unnest*, *cvt:pre:sort*).
  - Option *cvt:pre:boundlogarg* (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


## 20250801
- Changes in MP
  - Tolerances set by options *pre:feastol*,
    *pre:feastolrel* both need to be violated
    to produce a warning on contradicting
    variable/constraint bounds. Previously
    the preprocessor failed on any violation,
    without letting the solver try.
  - Options *cvt:compl*, *cvt:compl:eps* control
    complementarity reformulations.
  - Multi-objective emulator: added support for
    objective-specific options via objective suffixes
    beginning with *option_*
  - Option *cvt:unnest*: bits 2 and 4 switch on
    inlining of linear and quadratic subexpressions
    produced during reformulations (by default on).
  - Options *cvt:pre:ctx2ineq*, *cvt:pre:ctx2count*
    to control context propagation into conditional
    comparisons #267.


## 20250521

- Initial version
