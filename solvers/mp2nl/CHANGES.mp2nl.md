Summary of recent updates to MP2NL for AMPL
==========================================


## 20261029
- Configuration *alglib4nl* (ALGLIB for AMPL), from Knitro's: logical
  expressions, if-then-else, min/max/abs, SOS and indicators are
  linearized by MP; complementarity is converted with
  *cvt:compl=1* (product): ALGLIB ignores NL complementarity
  constraints, and the disjunctive and Fischer-Burmeister forms failed
  on the tests. Its time limit and relative gap are *timeout* and
  *pdgap*.
- *lim:time* (did nothing), and new *mip:gap*, *mip:gapabs*, are
  passed to the subsolver as its own keywords (e.g., Knitro's
  *maxtime*, *mip_opt_gap_rel*, *mip_opt_gap_abs*; BARON's *maxtime*,
  *epsr*, *epsa*) for the known configurations. They can be set per
  objective in a multi-objective solve (suffixes *.option_timelimit*,
  *.option_mipgap*), which used to fail with "Option not found".
- Fix: in a multi-objective solve (*obj:multi*), the NL file of
  every pass after the first had a wrong header (counts only for
  the items added since the previous pass); the subsolver read
  wrong data or stopped.
- Objective suffixes (e.g., *.objpriority*) are no longer passed
  on to the auxiliary variables and constraints of the objective's
  conversion. Knitro 16.0.0 with *use_asl=0* misreads binary NL
  files with any suffix, so a multi-objective solve of a model with
  logical expressions returned a zero solution.


## 20260624
- Changes in MP, in particular:
  - Option *cvt:pow2_as_qp* prefers quadratization
    of the ^2 operator for affine arguments (default).
  - Faster parsing of quadratic expressions with
    defined variables.
  - The main driver code is separated into a shared
    library and used commonly with AMPLS API.


## 20260525
- Forward AMPL imported function calls to the subsolver,
  such as gsl_ran_gaussian_pdf().


## 20260515
- Randomly permute algebraic constraints, which
  seems to help some solvers.
- Produce defined variables for algebraic expressions
  with reference count above *nl:assign:defvar*.
- Changes in MP, in particular:
  - Option *cvt:expr:nlreif*: above which reference count,
    logical expressions are reified.

## 20260414
- Changes in MP, in particular:
  - Option *cvt:pre:boundsbest* defaults to making
    auxiliary subexpression result variables free,
    unless the result is constant. This allows
    stronger nonlinear presolve.


## 20260108
- Reuse solutions between iterations of
  multi-objective emulator.


## 20251210
- Change in MP: option *obj:multi:options*: fix integer-valued
  objective-specific options.


## 20251121
- New option *obj:multi:options* to control
  whether multi-objective option suffixes
  are used.


## 20251021
- Changes in MP: option *cvt:expr:nlassign*.


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


## 20250128
- Updates for compilers compatibility


## 20241219
- Added *WantLogicalizedProd2Bin()*


## 20240926
- First release of the MP2NL driver
