Summary of recent updates to MP2NL for AMPL
==========================================


## unreleased
- Windows: a temporary directory with a long path (close to or beyond
  the 248-character limit for directories) is now supported. When the
  temporary directory cannot be used, MP2NL stops with an error message
  instead of terminating abnormally.

## 20261029
- Version 0.2.
- Native NL complementarity remains the default (*acc:nlcompl=2*),
  including with Knitro's configuration.
- Use *acc:nlcompl* to control native complementarity in MP2NL;
  *acc:compl* has no effect. Reformulations selected by *cvt:compl*
  apply only with *acc:nlcompl=0*.
- New configurations of solvers not accepting NL complementarity
  convert it (*acc:nlcompl=0*): *gurobiasl*, *cplexasl*, *xpressasl*
  (disjunction), *bonmin*, *snopt*, *minos*, *lgo*, *loqo* (product,
  *cvt:compl=1*), in addition to *ipopt*, *lindoglobal* and others.
  MINOS and SNOPT ignored NL complementarity
  constraints, the others failed.
- Quadratic constraints created by reformulations (e.g., *acc:div=0*)
  continue to be written as NL quadratic expressions by default.
  *cvt:quadcon=0* now requests linearization (previously ignored for NL
  output). ASL MIP solvers may need their nonconvex options
  (e.g., Gurobi's *nonconvex=2*).
- Configuration *alglib4nl* (ALGLIB for AMPL): logical
  expressions, if-then-else, min/max/abs, SOS and indicators are
  linearized by MP; complementarity is converted with
  *cvt:compl=1* (product): ALGLIB ignores NL complementarity
  constraints. The product form is the supported default for this
  configuration. Its time limit and relative gap are *timeout* and
  *pdgap*.
- *lim:time* (did nothing), and new *mip:gap*, *mip:gapabs*, are
  passed to the subsolver as its own keywords (e.g., Knitro's
  *maxtime*, *mip_opt_gap_rel*, *mip_opt_gap_abs*; BARON's *maxtime*,
  *epsr*, *epsa*) for the known configurations. They can be set per
  objective in a multi-objective solve (suffixes *.option_timelimit*,
  *.option_mipgap*), which used to fail with "Option not found".
- Fix: multi-objective passes after the first could submit an invalid
  NL model, causing the subsolver to read wrong data or stop.
- Fix: multi-objective models with logical expressions could return
  an incorrect zero solution with Knitro 16.0.0 and *use_asl=0*.
  Unnecessary objective suffixes are now omitted from submitted models;
  Knitro's binary NL reader still has limitations with suffixes.


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
