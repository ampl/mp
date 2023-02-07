Summary of recent updates to COPT for AMPL
==========================================


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
    

### 20220411
- First mp-based release, linked with COPT libraries 4.0.5
