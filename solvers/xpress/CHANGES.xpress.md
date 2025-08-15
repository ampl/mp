Summary of recent updates to Xpress for AMPL
============================================


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.
- Fixed modeling interface for AND, OR, MIN, MAX.
- Option alg:numericalemphasis.
- Handling Ctrl-C correctly.


## 20250806
- Updated to Xpress 45.01.02 (9.7), which includes
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


## 20250531
- Updated to Xpress 45.01.01 (9.6) that includes:
  - Pre root parallel heuristics phase.
  - New deterministic measure of algorithmic work.
  - Performance improvements.
- Option *lp:opttol* (*opttol*), synonym for previous
  *lp:optimalitytolerance*.
- Added keywords *lim:work*, *lim:prerootwork* and *mip:prerooteffort*.
- Added value *3* to *pre:domcol*.


## 20250429
- Fix a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250426
- Relinked with XPRESS 44.01.04.
- Non-linear constraints are now passed via the expression
  API, that can improve performance.
- Changes in MP
  - Option *cvt:qp2pass* (default even faster parsing
    of quadratics)


## 20250329
- Changes in MP:
  - Option *cvt:multoutcard* to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.
- Renamed option poolnbest to sol:poollimit.
	Old name kept as synonym.


## 20250308
- Fix several option descriptions.
- Default to natively accept nonlinear expressions
  since they seem to perform well in v9.5.
  Full XPRESS license required.


## 20241227
- Fix acceptance options for nonlinear formulas
  (still disabled by default.)
- Option *pre:solve_nlp*
- Updated to Xpress 9.5 (44.01.01):
  - Performance improvements in the Global optimization engine
    and in the MIP solver
  - New column-based presolve eliminations to reduce model size when free variables are present in the model instance.
  - In models with big-M coefficients, additional tightening of coefficients is performed.  
  - Knapsack cover inequalities are now strengthened by using superadditive lifting.


## 20240823
- Updated to Xpress 43.01.03
- Added hybrid gradient algorithm (set *bar:alg* to 4)
- Option *bar:cpuplatform* defaults to -2. To replicate the previous
  versions behaviour, set to -1.


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


## 20240606
- Fix mip:basis (post-solving as fixed MIP),
	although fix + resolve from AMPL can be faster.
- Fix reporting infeasibility with obj:multi=1.


## 20240604
- Presolve division by constant, resulting in fewer constraints
- Fix no-solution case in multi-objective emulator


## 20240529
- *Multi-objective emulator*
	- All flat MP solvers support multi-objective mode (obj:multi=1),
		either natively, or via emulation.
	- Suffixes .objpriority, .objweight, .objabstol, .objreltol.
	- [BREAKING] Default intuitive handling of .objweight,
		see option obj:multi:weight, even when natively supported.


## 20240429
- [BREAKING] Merged `report_times` and `timing`; they 
  are now aliases, set the value to 1 to have basic info,
  to 2 to have more detailed info.

  
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
    are passed natively to the solvers accepting them,
    vs previously quadratic or linear reformulation.
- *Option `report_times`*.
- *Unused `acc:` options*.
  - The constraint acceptance options `acc:...`
    for non-handled constraints are ignored
    (previously triggered error.)


## 20240312
- Updated Xpress to 42.01.5
- *Xpress Global Solver*
  - Native implementation of most supported non-linear
    functions. Not enabled by default.


## 20240306
- *Xpress Global Solver*
  - Apply global (MI)NLP solver for non-convex
    non-linear models.


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


## 20231117
- MP update: fixed graceful exit on Ctrl-C from AMPL in Linux
  and fixed issue with reading text-format NL files


## 20231103
- Improved translation of logical constraints:
  inlining of nested disjunctions and conjunctions;
  fewer auxiliary binary variables.


## 20231017
- Fixed a bug in NL reader on Windows.


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


## 20230808
- Updates for Xpress 9.2 (42.01.01), which include many performance
  improvements in the MIP solver. 
- Non-convex quadratic problems are now solved to global optimality 
  by default.
- Added options 'mip:heurshiftprop', 'tech:backgroundthreads' and
 'tech:backgroundselect'.
- The range of the control 'mip:feasibilityjump' was extended.


## 20230728
- Option 'tech:writesolution'.


## 20230727
- Fixed a memory leak #217
- Fixed quadratic constraints


## 20230726
- MP: fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230724
- Option [solver_]auxfiles rc; transfers names
	of variables and linear constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230714
- Options barrier/primal/dual/network like in ASL.


## 20230621
- Fix quadratic objective with repeated subexpressions.
- Fix reformulation of ==> / else.


## 20230616
- Changes in MP.
- Eliminated spurious warnings
- Fixed passing of quadratic objective


## 20230607
- Amended detection and display of licensing errors
- Licensing allows now usage of an XPRESS license from an AMPL-based
  deployment


## 20230603
- Added option 'tech:logfile' to enable output to a log file


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution

  
## 20230227
- Eliminated warning message when a non feasible solution is added as a starting
  point for the MIP search


## 20230207
- *Changes in the MP library*


## 20221228
- Changes in MP


## 20221222
- Bug fixes in MP


## 20221211
- *Changes in MP: added the ==> else operator*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *Changes in MP: PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.


## 20221208
- First release of the MP-based driver for Xpress (solver version 9.0)
