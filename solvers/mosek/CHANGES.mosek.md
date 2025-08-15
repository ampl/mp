Summary of recent updates to MOSEK for AMPL
===========================================


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


## 20250429
- Fix a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250426
- Updated to MOSEK 11.0.18
- Return solution whenever available (even 'unknown'
  or 'infeasible'). This seems helpful on some
  numerically tough models.
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.
- Option mip:conic:outapprox, useful to find
  feasible solutions in numerically tough conic
  models.


## 20250308
- Options mip:gapabs, mip:feastol, mip:heurlevel, mip:feaspump


## 20250226
- Updated to MOSEK 11.0.8.
- Option *bar:basis* to control crossover after interior-point
  optimization.
- Option *pre:scale* to control if to apply scaling during presolve
- On MacOS, Intel processors support is dropped

- 
## 20240901
- Option *tech:logfile* to write Mosek's output to a log file


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
- Updated to MOSEK 10.2, which provides many bugfixes.
- Added options `mip:varselection`, `pre:dualray_analysis` and 
  `lim:sol`.


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


## 20240118
- *Native Mosek options*
  - Use tech:optionnative, tech:optionnative(read/write)
    to read/write native options inline or to/from files.
    See mp.ampl.com/features-guide.html#solver-options.


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
- Improved translation of *SOCP vs QCP constraints*.
  - For convex QCP problems, don't attempt SOCP form.
    Previously some QCP problems required cvt:socp=0.
  - For conic problems, substitute QP terms
    from the objective as an auxiliary SOCP constraint
    (only for simple QP terms, otherwise should be
    done manually.)
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
- Print a warning when not all quadratic constraints
    have been converted to cones.


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
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Fixed a bug where equivalent conditional
  comparisons were not unified.


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230724
- option [solver_]auxfiles rc; transfers names
	of variables and linear constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230621
- Fix quadratic objective with repeated subexpressions.


## 20230616
- Changes in MP.


## 20230531
- Cones: recognize (affine_expr) >= y * exp(z/y)
  as exponential cone.
- Cones: recognize xy >= 1 as rotated SOC.
- Wrong solver options are gracefully reported via
  solve_message.


## 20230515
- *Exponential cones*. MP driver recognizes exponential
  cones from their algebraic representation and passes
  them to Mosek.


## 20230505
- *Updated Mosek library* to version 10.0.43. It includes a 
  number of bug fixes, including a numerical problem that 
  could occur with disjunctive constraints


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution
  

## 20230330
- Removed the message on missing dual solution for MIP.


## 20230321
- Recognition of second-order conic constraints
  from algebraic representations.

- First release in the AMPL bundle.


## 20220420
- First release of mock driver.
