Summary of recent updates to COPT for AMPL
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
- Enable native indicator constraints by default
  (options acc:ind..)


## 20250429
- Fix a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250424
- Changes in MP
  - Option cvt:qp2pass (default even faster parsing
    of quadratics)


## 20250408
- Fixed a bug in parsing of quadratic expressions


## 20250407
- Updated to COPT version 7.2.6
- GPU is now supported for PDLP on Windows and Linux;
  CUDA library version >= 11.7 is required. See
  options *lp:method* and *lp:pdlpgpumode*.


## 20250329
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.


## 20250308
- Option mip:gapabs


## 20241228
- Updated to COPT 7.2.4, which includes performance improvements to 
  the MIP, SOCP and QCQP solvers


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
- Updated to Copt 7.1.3, which provides many bugfixes.


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


## 20240410
- Updated to Copt 7.1.1, with support for first-order 
  method (PDLP). On Windows this can be used with GPU
  supporting CUDA.
  It can be controlled setting `lp:method` to the new 
  value `6` (PDLP) and with options `lp:pdlpgpumode`, 
  `lp:pdlpgpudevice` and `lp:pdlptol`


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


## 20240115
- *Solve result codes*.
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
- Updated to Copt 7.0.3, which includes many performance improvements
- New keyword `lim:soltime` to specify a time limit after a solution
  has been found
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


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230724
- Option [solver_]auxfiles rc; transfers names
	of variables and linear constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230621
- Fix quadratic objective with repeated subexpressions.


## 20230616
- Changes in MP.


## 20230428
- Updated to Copt 6.5.2, which includes performance improvements and
  fixes for issues with MIP problems


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution
  

## 20230330
- Recognition of second-order conic constraints
  from algebraic representations and conversion into
  quadratic constraints; COPT appears to recognize
  second-order cones from quadratics.


## 20230207
- *Changes in the MP library*
- Updated to Copt 6.0.4, which includes bugfixes
- Added support for unbounded/Farkas rays calculation


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


## 20230207
- Updated to Copt 6.0.1, which includes:
     many performance improvements 
     native support for quadratic constraints


## 20221012
- *Piecewise-linear approximation of quadratics*
    For non-convex quadratics, set the following options:
    cvt:quadobj=0 and/or cvt:quadcon=0.


## 20220928
- *Changes in MP*: piecewise-linear approximations of nonlinear functions,
    default value of big-M


## 20220715
- Updated to Copt 5.0.1, which includes many performance improvements
- Added feasibility relaxation (see *alg:feasrelax*)
- New parameters: *alg:iismethod*


## 20220615
- New parameter: *crossover*
- Minor changes to parmeter names


## 20220526
- *SOS constraints* are now detected also if the .ref suffix is integer
- Minor changes to parmeter names

## 20220511
- *Complementarity constraints: also quadratics*
    Complementarity constraints now handle quadratics.

- *Branch develop is used for new code*
    The active development branch is now *develop*.

- *Convert quadratic range constraints to QuadCon(LE/EQ/GE)*
    COPT does not support quadratic range constraints.
    Conversion of linear range constraints into one-side rhs
    constraints has been generalized for any algebraic ones.
    

## 20220411
- First mp-based release, linked with COPT libraries 4.0.5
