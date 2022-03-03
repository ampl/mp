Summary of recent updates to x-gurobi for AMPL
==============================================

### 20220303
- *Fix strict comparison tolerance*
    Option *cvt:mip:eps*, default 1e-3.

### 20220217
- *See also MP Library changes 20220216*

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
