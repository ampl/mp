Summary of recent updates to CBCMP for AMPL
===========================================

## 20230817
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Fixed a bug where equivalent conditional
  comparisons were not unified.


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230713
- Fixed how certain options are passed to the underlying CBC library

- Fixed returning of constraint duals


## 20230616
- Changes in MP.


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution

  
## 20230207
- *Changes in the MP library*


## 20221228
- *Changes in the MP library*


## 20221222
- *Changes in the MP library*


## 20221211
- *Changes in MP: added the ==> else operator*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *Changes in MP: PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.


## 20221208
- First release of MP-based cbc driver
