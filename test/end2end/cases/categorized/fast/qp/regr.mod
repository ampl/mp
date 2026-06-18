##########################
# Origin: sp #84

# option presolve 0;
option randseed 42;

param NA default 30;
param NS default 50;

set ATTRIBUTES := 1..NA;
set SAMPLES := 1..NS;

param mmX default 0.05;
param vmX default 0.0025;
param mvX default 0.05;
param vvX default 0.0004;

param mW default 1/NA;
param vW default 0.16;

param mE default 0.13;
param vE default 13;

param mX {ATTRIBUTES} default Normal(mmX, vmX);
param vX {ATTRIBUTES} default Normal(mvX, vvX);

param W {ATTRIBUTES} default Normal(mW, vW);
param X {a in ATTRIBUTES, SAMPLES} default Normal(mX[a], vX[a]);
param Y {s in SAMPLES}
   default
      sum {a in ATTRIBUTES} X[a,s] * W[a]
      + Normal(mE,vE);

var weight {ATTRIBUTES} >= 0;

var error_sample {d in SAMPLES} =
   Y[d] - sum {b in ATTRIBUTES} (weight[b] * X[b, d]);
var regularizer =
   sum {b in ATTRIBUTES}
      (weight[b] - 1 / card(ATTRIBUTES)) ^ 2 / 900;

minimize TotalPenalty:
    sum {d in SAMPLES} error_sample[d]^2 + regularizer;

subject to Normalize:
   sum {b in ATTRIBUTES} weight[b] = 1;