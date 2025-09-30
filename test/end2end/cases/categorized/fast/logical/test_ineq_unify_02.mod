
###########################################
## Test unification of (in)equalities
## test_ineq_unify_02.mod: on real-valued expressions.
## b[i] >< 0.5 converted to just b[i] or !b[i]
###########################################

param cmpEps default 1e-4;  # cvt:mip:eps

var b {1..19} binary;

var x {1..5} >=0 <=5;

minimize Obj:
  sum {i in 1..19} b[i]
    - x[2] - 2*x[4] + x[1] - 5*x[3] - x[5];

## ineq2related
s.t. C1: b[1]>0.5 ==> x[2] - x[4] >= -2;
s.t. C2: b[2]<0.5 ==> x[2] - x[4] < -2;
s.t. C3: x[2] - x[4] > -2-cmpEps || x[5] >= 3;
s.t. C4: x[2] - x[4] <= -2-cmpEps ==> x[5] > 3-cmpEps;

# negated
s.t. C1_1: b[5]<0.5 ==> -x[2] + x[4] <= 2;

## ineq2bndeq
s.t. CB1: b[3]>0.5 ==> x[3] - 3*x[5] >= -15;
s.t. CB2: b[10]>0.5 ==> x[3] - 3*x[5] <= -15;
s.t. CB3: b[4]>0.5 ==> x[3] - 3*x[5] <= -15 + cmpEps/2;
s.t. CB4: b[5]>0.5 ==> x[3] - 3*x[5] > -15 + cmpEps/3;
                       ## Also < -15?
s.t. CB5: b[6]>0.5 ==> x[3] - 3*x[5] < -15 + cmpEps;
s.t. CB6: b[7]>0.5 ==> x[3] - 3*x[5] >= -15 + cmpEps;
s.t. CB7: b[7]<0.5 ==> x[3] - 3*x[5] >= -15 + cmpEps/2;
s.t. CB8: b[8]>0.5 ==> x[3] - 3*x[5] == -15;
s.t. CB9: b[9]>0.5 ==> x[3] - 3*x[5] != -15;

s.t. CB10: b[11]>0.5 ==> x[3] - 3*x[5] >= 5;
s.t. CB11: b[12]<0.5 ==> x[3] - 3*x[5] <= 5;
s.t. CB12: b[13]>0.5 ==> x[3] - 3*x[5] >= 5 - cmpEps/2;
s.t. CB13: b[14]>0.5 ==> x[3] - 3*x[5] < 5 - cmpEps/3;
                       ## Also > 5?
s.t. CB14: b[15]>0.5 ==> x[3] - 3*x[5] > 5 - cmpEps;
s.t. CB15: b[16]>0.5 ==> x[3] - 3*x[5] <= 5 - cmpEps;
s.t. CB16: b[17]<0.5 ==> x[3] - 3*x[5] <= 5 - cmpEps/2;
s.t. CB17: b[18]>0.5 ==> x[3] - 3*x[5] == 5;
s.t. CB18: b[19]>0.5 ==> x[3] - 3*x[5] != 5;
