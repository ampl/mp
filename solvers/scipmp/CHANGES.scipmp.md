Summary of recent updates to SCIP for AMPL
==========================================

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
