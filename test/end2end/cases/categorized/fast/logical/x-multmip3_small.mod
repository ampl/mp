## Trying to replicate the cvt:cmp:eps problem from x-multmip3.mod.
## Checks that all cases of x[i] >= minl are assigned
## the same auxiliary variable, otherwise with cvt:cmp:eps=0
## the optimal objective becomes 0. See also #102.

param n default 4;
param minl default 375;
param limU default 500;
param fcost default 50;
param overall default 1200;
param maxserve default n-1;

var x{1..n} >= 0;

minimize Total:
   sum {i in 1..n} if x[i]>=minl then fcost;

s.t. Lin01:
   sum {i in 1..n} x[i] >= overall;

s.t. Disj{i in 1..n}:
   x[i]==0  or minl <= x[i] <= limU;

s.t. Count:
   count {i in 1..n} (x[i] >= minl) <= maxserve;
