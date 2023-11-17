Summary of recent updates to COPT for AMPL
==========================================

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
