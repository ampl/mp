# -------------------------------------------------------------
# IIS LogicalCon 06
# infeas_lcon_06.mod
# Test that AMPL can receive suffixes when all cons logical ...
# AND there is no objective but a linear constraint.
# Test also that AMPL sends suffixes then.
# Modified example from solvers-private/#63.
# Addresses escrow/#141.
# -------------------------------------------------------------

set TYPE;
param T;
set TIME_STEP := 1..T ordered;
var elec_flow{TYPE, TIME_STEP} >= 0, <= 100000;
var power{TYPE} >= 0, <= 100000;
param min_load{TYPE} >= 0;

s.t. flow {i in TYPE, j in TIME_STEP}:
    elec_flow[i,j] > 0
      ==> 15 + min_load[i]*power[i] <= elec_flow[i,j] <= power[i];

s.t. LC1: sum {i in TYPE} power[i] <= 5000;

data;
param T := 2;
set TYPE := FC ELLY;
param min_load := FC 1.2 ELLY 2;

model;
drop flow ['ELLY', 2];

s.t. Flow4 {i in {'ELLY'}, j in {2}}:
    (elec_flow['FC', 1] == 0
    && elec_flow['ELLY', 1] == 0
    && elec_flow['FC', 2] == 0)
    ==> 16 + min_load[i]*power[i] <= elec_flow[i,j] <= power[i] + log(power[i]);

suffix funcpieces IN;

let Flow4['ELLY', 2].funcpieces := -2;
let flow['ELLY', 1].funcpieceratio := 1e-3;
