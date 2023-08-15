Summary of recent updates to SCIP for AMPL
==========================================

## 20230815
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Changed default tolerance for strict comparisons
  to 0 (option cmp:eps, #102.)
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
