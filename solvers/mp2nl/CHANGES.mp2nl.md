Summary of recent updates to MP2NL for AMPL
==========================================


## unreleased

- The result line now starts with MP2NL's name and then the subsolver's,
  for example "MP2NL 0.2: Gurobi 14.0.0: optimal solution", each once.
  Before, the subsolver's name was printed before MP2NL's, so that the line
  read "Gurobi 14.0.0: MP2NL 0.2: Gurobi 14.0.0: optimal solution" (with
  ASL-based subsolvers; with MP-based ones, wherever the output is not a
  terminal).
- The debug statistics returned by the subsolver (problem suffixes
  `flat0_*`, `flat1_*` and `option___*` with `debug=1`) are now reported
  with prefix `sub_`, for example `sub_flat1__lineq`. Before, they
  overwrote MP2NL's own statistics of the same names.
- Fix: a solution file left by an earlier solve (a reused `nl:stub`, or the
  previous pass of a multi-objective emulation) was read as the result of a
  solve in which the subsolver wrote no solution. Such a solve now fails with
  a message. If the old file cannot be removed, the driver reports an error
  before launching the subsolver.
- Fix: a solution written by the subsolver before it exited with a nonzero
  code was discarded; it is kept now, and the exit code is reported in a
  warning.
- Fix: subsolver options (`solver_options`) were interpreted by the shell:
  characters such as `;` or `$` acted as commands and quoted values were
  split. Each option is now passed as one argument; quotes group words.
  An unclosed quote is an error.
- Fix: the solution values, duals and suffixes of a previous pass were
  reported after a pass that failed or returned no such data.
- Fix: logical constraints were counted as nonlinear objectives in the NL
  header.
- Fix: NL files written through the NLModel API for quadratic objectives with
  one or two terms were invalid, and nonlinear variables were miscounted
  (CONOPT returned a wrong optimum).
- When the temporary directory cannot be used (missing, or on Windows
  a path close to the 248-character limit for directories), MP2NL stops
  with an error message instead of terminating abnormally:
  the failure and its message are reported to AMPL in the solution file.
- The NL file stub is quoted on the solver command line, so temporary
  directories with spaces or other shell characters in their path work.
  Exception on Windows: '%' is still expanded by the command interpreter,
  so a path containing a variable reference such as `%OS%` fails.
- Library API: the close and model functions are now *AMPLSClose_mp2nl*
  and *AMPLSGetModel_mp2nl*, lower case like *AMPLSOpen_mp2nl* and the
  other MP libraries (were *AMPLSClose_MP2NL*, *AMPLSGetModel_MP2NL*).
- Fix: a subsolver solution with fewer duals than algebraic
  constraints could corrupt memory.


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
