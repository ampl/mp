Summary of recent updates to MOSEK for AMPL
===========================================

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
