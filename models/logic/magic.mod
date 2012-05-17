
param n integer > 0;

var X {0..n} integer >= 0, <= n - 1;

subj to magicality {j in 0..n}:
   X[j] == numberof j in ({k in 0..n} X[k]);

subj to redundant:
   sum {j in 0..n} j * X[j] = n + 1;
