Summary of recent updates to SCIP for AMPL
==========================================


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


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
- Reenable native OR constraint by default
  (requires marking all variables as integer,
  even if fixed). See option acc:or.
- Enable native indicator constraints and cos().
  Options acc:ind.. and acc:cos.


## 20250429
- Fix a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250426
- Relinked with SCIP version 9.2.2
- Disabled native OR constraint by default
  - Enabled native AND constraint which seems
    correct now


## 20250426
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.


## 20250404
- Fixed:
  - Option lp:threads
  - Option lp:presolving
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.
- Option sol:poollimit.


## 20250308
- Changes in MP.


## 20250204
- Updated to SCIP 9.2
- MINLP expression trees (option acc:_expr.)
	- Using the expression tree API, available since SCIP v8,
		to model nonlinear expressions. Earlier the expressions
		were submitted to the solver in a flat form equating
		an auxiliary variable to each expression.


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
- Presolve division by constant, resulting in fewer constraints
- Fix no-solution case in multi-objective emulator


## 20240531
- Updated to SCIP 9.0.1


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
    pl_linearize 1`, default; set to 0 to use 
    MP linearization.)
  - Disallow repeated weights for SOS constraints
    (suffixes `.sosno`/`.ref`.)
- *Option `report_times`* 
- *Unused `acc:` options*.
  - The constraint acceptance options `acc:...`
    for non-handled constraints are ignored
    (previously triggered error.)


## 20240221
- Updated to SCIP 9.0


## 20240121
- Updated to SCIP 8.1.0


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


## 20231015
- Added aliases mip:gap and mip:gapabs for consistency with
  other solvers


## 20231012
- Fixed bug in option tech:logfile (logfile).


## 20239021
- Updated to SCIP 8.0.4 and SoPlex 6.0.4


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


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230625
- Fix dual solutions.


## 20230623
- Added support for quadratic cone constraints but not recommended.


## 20230622
- Added NumQPCons, NumSOSCons and NumIndicatorCons functions.


## 20230621
- Fix quadratic objective with repeated subexpressions.


## 20230619
- Added MULTISOL support.


## 20230616
- Changes in MP.


## 20230515
- First release of mock driver
