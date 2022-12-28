Summary of recent updates to gurobi for AMPL
==============================================

## 20221228
- Changes in MP


## 20221222
- *Fixed #195": shorter error message for missing NonConvex=2 option

- Bug fixes in MP


## 20221211
- *Changes in MP: added the ==> else operator*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *Changes in MP: PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.


## 20221113
- *Gurobi 10.0 support*
    Options 'lim:mem', 'lim:softmem', 'mip:obbt' and  'alg:networkalg'
    Added option tech:writepresolvedprob to export the presolved model
- *Released the new MP-based Gurobi driver*
    The new driver becomes the default (and is named just 'gurobi')

## 20221012
- *Piecewise-linear approximation of quadratics*
    For Gurobi, non-default.
    To use, set the options cvt:quadobj=0 cvt:quadcon=0.
    Recognizing x^2 for stronger univariate approximation

## 20220928
- *Changes in MP*: piecewise-linear approximations of nonlinear functions,
    default value of big-M

- For *range constraints* x-gurobi reports nonbasic status low/upp,
    for one-sided constraints low/upp/equ, consistent with ASL drivers

- *Sensitivity analysis*: use constraint suffixes .sens(lb/ub)(lo/hi),
    the old-style suffixes .sensrhs(lo/hi) meaningful only for one-sided constraints.

- Suffixes *.iis(lb/ub)force* on constraints and variables

## 20220802 
- Added support for 'params' option from command line and environment variable

## 20220725
- *Changes in MP* fixed suffixes export on Windows and multiple solutions handling

## 20220720
- *Options 'funcpieces', 'funcpiecelength', 'funcpieceratio', 'funcpieceerror'*
    The above options (and corresponding suffixes) are passed to Gurobi. The
    suffixes can specialize the values for individual constraints.

    Subexpressions: note that if a subexpression is contained in several
    constraints, for contradicting suffix values the maximum is taken.


## 20220706
- *Relinked with Gurobi 9.5.2, which contains bug fixes*


## 20220511
- *Complementarity constraints: also quadratics*
    Complementarity constraints now handle quadratics.

- *Branch develop is used for new code*
    The active development branch is now *develop*.

- *Convert quadratic range constraints to QuadCon(LE/EQ/GE)*
    Gurobi does not support quadratic range constraints.
    Conversion of linear range constraints into one-side rhs
    constraints has been generalized for any algebraic ones.


### 20220408
- *Linear complementarity in MP: 1st go*

- *Other changes in MP*


### 20220303
- *Fix strict comparison tolerance*
    Option *cvt:mip:eps*, default 1e-3.


### 20220217
- *Assume new constraints are active (#152)*:
    Gurobi requires a complete basis for hotstart and we have to guess the
    statuses of new variables and constraints.
    
- *Allowing SOS constraints with repeated weights (#163)*:
    Although Gurobi states SOS weights should be unique, it accepts them repeated.
    This happens when AMPL linearizes a PL function with redundant (repeated) slopes.
    It seems better to use PL functions natively (*option pl_linearize 0;*).
    
- *Native handling of abs, min/max, and/or, and indicators by default*:
    For the general constraints abs, min/max, and/or, Gurobi 9.5 seems to use
    tight MIP reformulations, matching the performance of MIPConverter redfinitions.
    In contrast, indicator constraints behave differently to MIP reformulations
    (accessible by acc:ind_..=1): better primal and worse dual bounds.
    Setting acc:* = 2 as default (native handling).


### 20220202
- *Basis status low/upp/sup for new variables*:
    when new variables are added, AMPL assigns .sstatus *none* while Gurobi 9.5 
    needs a complete basis so we automatically set Gurobi var status to *low*/*upp*/*sup*
    depending on where 0.0 is relative to the bounds.


### 20220128
- First eXperimental release, linked with Gurobi 9.5.
