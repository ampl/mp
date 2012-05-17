
param n integer > 0;
var Row {1..n} integer >= 1 <= n;

subj to c1: alldiff ({j in 1..n} Row[j]);
subj to c2: alldiff ({j in 1..n} Row[j]+j);
subj to c3: alldiff ({j in 1..n} Row[j]-j);
