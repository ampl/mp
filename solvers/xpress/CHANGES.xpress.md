Summary of recent updates to Xpress for AMPL
============================================

## 20230726
- Fixed inequalities of integer expressions with
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
