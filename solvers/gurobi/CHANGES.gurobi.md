Summary of recent updates to gurobi for AMPL
============================================


## 20250617
- Option mip:opttol renamed as lp:opttol (opttol).
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


## 20250515
- Updated to Gurobi 12.0.2, which includes bugfixes
- Improved detection of unsupported multiobjective models


## 20250429
- Added the Tan expression (previously always submitted
	to Gurobi as a general constraint).
- Fixed a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250424
- Changes in MP
  - Option cvt:qp2pass (default even faster parsing
    of quadratics)


## 20250329
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.


## 20250311
- Additional diagnostic messages are printed in case the solver
  fails to initialize


## 20250308
- Option alg:kappa_exact


## 20250204
- Updated to Gurobi 12.0.1


## 20241119
- MP changes: added support for x^y expressions, where x and y
  are variables


## 20241114
- Updated to Gurobi 12.0 which provides performance improvements 
  across a variety of model families
- Gurobi now supports nonlinear constraints of the form y=f(x),
  where f is a multivariate function given by a closed-form expression.
- Solve mixed-integer nonlinear programming (MINLP) problems 
  to global optimality more efficiently via the new non-linear 
  expression interface.


## 20240808
- Added option *alg:nlpheur* to control an heuristic for
  non-convex quadratic models.


## 20240728
- Updated to Gurobi libraries 11.0.3, which include many  
  bug fixes.

 
## 20240724
- Option *acc:_all*
	- Useful to disable all reformulations (acc:_all=2),
		or force linearization (acc:_all=0).
- Option *cvt:prod*     
  - Controls reformulation of binary products into logical 
    constraints.
- Faster input of quadratic expressions.


## 20240617
- *Multi-objective emulator*
	- obj:multi=2 forces emulation, even if MO natively supported.
	- Fixed a bug in the objective degradation suffixes
		.objasbtol, .objreltol.


## 20240604
- Presolve division by constant, resulting in fewer constraints.
- Fix no-solution case in multi-objective emulator.


## 20240529
- *Multi-objective emulator*
	- All flat MP solvers support multi-objective mode (obj:multi=1),
		either natively, or via emulation.
	- Suffixes .objpriority, .objweight, .objabstol, .objreltol.
	- [BREAKING] Default intuitive handling of .objweight,
		see option obj:multi:weight, even when natively supported.


## 20240518
- Updated to Gurobi libraries 11.0.2, which include many 
  bug fixes.


## 20240429
- [BREAKING] Merged `report_times` and `timing`; they 
  are now aliases, set the value to 1 to have basic info,
  to 2 to have more detailed info.

  
## 20240327
- Added support for Web License Server parameters, via options
  `wls_licenseid`, `wls_accessid`, `wls_secret`, `wls_token`, 
  `wls_tokenduration` and `wls_tokenrefresh`.


## 20240320
- *SOS constraints*.
  - Fixed handling of SOS2 constraints created by AMPL
    as reformulations of PL expressions (`option
    pl_linearize 1`, default; set to 0 to use the solver's
    native PL functions or MP linearization.)
  - Disallow repeated weights for SOS constraints
    (suffixes `.sosno`/`.ref`.)
- *Native handling of POW(x, INT)*.
  - Power expressions with positive integer exponent
    are passed natively vs previously quadratic or linear 
    reformulation. For best performance, global solving capability
    might be needed (`global=1`).
- *Option `report_times`* 
- *Unused `acc:` options*.
  - The constraint acceptance options `acc:...`
    for non-handled constraints are ignored
    (previously triggered error.)


## 20240311
- Added option `tech:reportwork` to display or return in the 
  problem suffix `work` the work units spent while solving 
  the problem.


## 20240310
- Updated to Gurobi libraries 11.0.1, which include many 
  bug fixes.


## 20240115
- *Solve result codes*
  - List codes by running (solver) -!
  - [BREAKING] Standardized codes. Major changes:
    - 100-199 (solved?) means solution candidate
      provided, but can be suboptimal/infeasible
    - 300-349 means unbounded problem but
      feasible solution returned
    - 400-449 means limit/interrupt but feasible
  - [BREAKING] sol:chk:fail returns code 150 (solved?)
- Improved translation of *SOCP constraints*.
  - Options cvt:socp, cvt:socp2qc.
- Compact solution check warnings
- Fixed presolve of the power function #226.


## 20231208
- Option `global` and suffix `.global` synonyms for
  `funcnonlinear`
- Changed values for `global` / `funcnonlinear` options
  and suffixes: default 0 (Gurobi automatic),
  -1 (static PL approximation), 1 (global solving).


## 20231206
- Updated to Gurobi 11
  - Non-linear models can now be solved using spatial 
    branch-and-bound and outer approximation instead of piecewise
    linearizattion. To use this set the option `pre:funcnonlinear` 
    to 1 and the suffix `funcnonlinear` to -1
  - Added keywords `cut:mixingcuts`,  `alg:concurrentmethod` 
    and `alg:solutiontarget`


## 20231117
- Added option lim:work (worklimit).
- MP update: fixed graceful exit on Ctrl-C from AMPL in Linux
  and fixed issue with reading text-format NL files


## 20231103
- Improved translation of logical constraints:
  inlining of nested disjunctions and conjunctions;
  fewer auxiliary binary variables.


## 20231017
- Fixed a bug in NL reader on Windows.


## 20230920
- Updated to Gurobi libraries 10.0.3, which include many bugfixes.


## 20230919
- *mp_options*.
	Receive mp_options from AMPL (for all MP solvers).
	They are parsed before (solvername)_options.
- Solution checking: relative tolerance
	sol:chk:feastolrel; options sol:chk:round, sol:chk:prec.


## 20230831
- Solution checking, options sol:chk:* (experimental).
- Preprocess And/Or constraints.


## 20230817
- Alternative solutions: solve status equal to that
  of the final solution.
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Fixed a bug where equivalent conditional
  comparisons were not unified.


## 20230728
- Option 'tech:writesolution' #218
- Option 'writeprob' ('tech:writemodel') ASL-compatible


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.
- Option 'writesol=filename' to export
  solutions/results.


## 20230724
- Option [solver_]auxfiles rc; transfers names
	of variables and constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230625
- Updated to Gurobi libraries 10.0.2, which include many bugfixes.


## 20230621
- Fix quadratic objective with repeated subexpressions.


## 20230616
- Changes in MP.


## 20230531
- Wrong solver options are gracefully reported via
  solve_message.


## 20230522
- Added option 'lim:sol' to set a limit on the number of solutions found


## 20230426
- Fixed partial MIP start.


## 20230330
- Fixed the outlev option.


## 20230321
- Recognition of second-order conic constraints
  from algebraic representations and conversion into
  quadratic constraints; Gurobi appears to recognize
  second-order cones from quadratics.


## 20230207
- *Changes in the MP library*.


## 20230206
- Relinked with Gurobi version 10.0.1, which includes many bugfixes

- Changed behaviour of 'tech:logfile', which no longer implies 'tech:outlev=1'.
  Specifying a logfile will enable gurobi's full log to file only; to obtain both 
  console and log output, set also 'tech:outlev=1'.


## 20221228
- Changes in MP

## 20221222
- *Fixed #195*: shorter error message for missing NonConvex=2 option

- Bug fixes in MP


## 20221211
- *Changes in MP: added the ==> else operator*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *Changes in MP: PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.


## 20221113
- *Gurobi 10.0 support*
    Options 'lim:mem', 'lim:softmem', 'mip:obbt' and  'alg:networkalg'
    Added option tech:writepresolvedprob to export the presolved model
- *Released the new MP-based Gurobi driver*
    The new driver becomes the default (and is named just 'gurobi')

## 20221012
- *Piecewise-linear approximation of quadratics*
    For Gurobi, non-default.
    To use, set the options cvt:quadobj=0 cvt:quadcon=0.
    Recognizing x^2 for stronger univariate approximation

## 20220928
- *Changes in MP*: piecewise-linear approximations of nonlinear functions,
    default value of big-M

- For *range constraints* x-gurobi reports nonbasic status low/upp,
    for one-sided constraints low/upp/equ, consistent with ASL drivers

- *Sensitivity analysis*: use constraint suffixes .sens(lb/ub)(lo/hi),
    the old-style suffixes .sensrhs(lo/hi) meaningful only for one-sided constraints.

- Suffixes *.iis(lb/ub)force* on constraints and variables

## 20220802 
- Added support for 'params' option from command line and environment variable

## 20220725
- *Changes in MP* fixed suffixes export on Windows and multiple solutions handling

## 20220720
- *Options 'funcpieces', 'funcpiecelength', 'funcpieceratio', 'funcpieceerror'*
    The above options (and corresponding suffixes) are passed to Gurobi. The
    suffixes can specialize the values for individual constraints.

    Subexpressions: note that if a subexpression is contained in several
    constraints, for contradicting suffix values the maximum is taken.


## 20220706
- *Relinked with Gurobi 9.5.2, which contains bug fixes*


## 20220511
- *Complementarity constraints: also quadratics*
    Complementarity constraints now handle quadratics.

- *Branch develop is used for new code*
    The active development branch is now *develop*.

- *Convert quadratic range constraints to QuadCon(LE/EQ/GE)*
    Gurobi does not support quadratic range constraints.
    Conversion of linear range constraints into one-side rhs
    constraints has been generalized for any algebraic ones.


## 20220408
- *Linear complementarity in MP: 1st go*

- *Other changes in MP*


## 20220303
- *Fix strict comparison tolerance*
    Option *cvt:mip:eps*, default 1e-3.


## 20220217
- *Assume new constraints are active (#152)*:
    Gurobi requires a complete basis for hotstart and we have to guess the
    statuses of new variables and constraints.
    
- *Allowing SOS constraints with repeated weights (#163)*:
    Although Gurobi states SOS weights should be unique, it accepts them repeated.
    This happens when AMPL linearizes a PL function with redundant (repeated) slopes.
    It seems better to use PL functions natively (*option pl_linearize 0;*).
    
- *Native handling of abs, min/max, and/or, and indicators by default*:
    For the general constraints abs, min/max, and/or, Gurobi 9.5 seems to use
    tight MIP reformulations, matching the performance of MIPConverter redfinitions.
    In contrast, indicator constraints behave differently to MIP reformulations
    (accessible by acc:ind_..=1): better primal and worse dual bounds.
    Setting acc:* = 2 as default (native handling).


## 20220202
- *Basis status low/upp/sup for new variables*:
    when new variables are added, AMPL assigns .sstatus *none* while Gurobi 9.5 
    needs a complete basis so we automatically set Gurobi var status to *low*/*upp*/*sup*
    depending on where 0.0 is relative to the bounds.


## 20220128
- First eXperimental release, linked with Gurobi 9.5.
