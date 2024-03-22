Summary of recent updates to HiGHS for AMPL
===========================================


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


## 20240307
- *Updated* to HiGHS 1.7.0
- Added solver pdlp (use option `alg:method` to specify it and
  `alg:pdlpdgaptol`, `lim:pdlpnativetermination`, `pre:pdlpscaling` 
  and `lim:pdlpiterationlimit`).
- Added optional pre centring steps (see options `pre:centring`,
  `pre:maxcentringsteps` and `pre:centringratiotolerance`).
- Added options `pre:userboundscale`, `pre:usercostscale`, 
  `lim:objectivebound` and `lim:objectivetarget`.


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
- Compact solution check warnings
- Fixed presolve of the power function #226.


## 20231117
- MP update: fixed graceful exit on Ctrl-C from AMPL in Linux
  and fixed issue with reading text-format NL files


## 20231109
- *Updated* to HiGHS 1.6.0

## 20231103
- Improved translation of logical constraints:
  inlining of nested disjunctions and conjunctions;
  fewer auxiliary binary variables.


## 20231017
- Fixed a bug in NL reader on Windows.
- Don't return basis for MIPs.
    Caused trouble with older versions of AMPL
    with logical constraints (workaround: alg:basis=0.)


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


## 20230728
- Option 'tech:writesolution'


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.


## 20230724
- option [solver_]auxfiles rc; transfers names
	of variables and linear constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230621
- Fix quadratic objective with repeated subexpressions.
- Fix Hessian API.


## 20230616
- Imitate partial MIP start.
- Changes in MP.


## 20230531
- *MIP start*.
  HiGHS 1.5 supports complete MIP starts only.
- Wrong solver options are gracefully reported via
  solve_message.


## 20230522
- Fixed solution status reporting to AMPL
- Fixed basis input with obj offset or missing data
- Added warm start for LP problems
- Reading column names from *col* file is present


## 20230424
- *Changes in the MP library*: added variable names support
  and removed spurious starting solution
  

## 20230227
- Fixed a problem when retrieving basis status
- Implemented retrieval of MIP gap so that if optimality is reached via presolving,
  the returned gap is 0


## 20230224
- *Updated* to HiGHS 1.5.1; now returning number of LP iterations in a MIP solve


## 20230222
- Fixed: now returning the correct dual values for the constraints


## 20230209
- *Updated* to HiGHS 1.4.2
- Fix for error encountered when passing variables status if variables with unknown 
  basis status are present.


## 20230207
- *Changes in the MP library*


## 20221228
- Changes in MP


## 20221222
- *Updates for HiGHS 1.4.1*
    - The *run_crossover* option has now values "on", "off" and "choose". The latter 
      results in crossover being run if the result of IPM without crossover is imprecise
    - Bug fixes

## 20221211
- *Changes in MP: added the ==> else operator*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *Changes in MP: PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.


## 20221012
- *Piecewise-linear approximation of quadratics*
    HiGHS accepts quadratic objectives.
    For nonconvex ones, set cvt:quadobj=0.
    Recognizing x^2 for stronger univariate approximation


## 20220928
- *Changes in MP*: piecewise-linear approximations of nonlinear functions,
    default value of big-M


## 20220603
- Fixed an issue with passing/retrieving basis information


## 20220524
- First release of HiGHS for AMPL
