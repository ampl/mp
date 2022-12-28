Summary of recent updates to HiGHS for AMPL
===========================================

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


### 20220603
- Fixed an issue with passing/retrieving basis information


### 20220524
- First release of HiGHS for AMPL
