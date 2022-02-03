Summary of recent updates to x-gurobi for AMPL
==============================================

### 20220202
- *Basis status low/upp/sup for new variables*:
    when new variables are added, AMPL assigns .sstatus *none* while Gurobi 9.5 
    needs a complete basis so we automatically set Gurobi var status to *low*/*upp*/*sup*
    depending on where 0.0 is relative to the bounds

### 20220128
- First eXperimental release, linked with Gurobi 9.5.
