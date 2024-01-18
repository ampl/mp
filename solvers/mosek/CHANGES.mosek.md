Summary of recent updates to MOSEK for AMPL
===========================================


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
