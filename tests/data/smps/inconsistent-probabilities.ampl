option auxfiles rc;
suffix stage IN;
set S = 0..1;
param P{S} = 1 / 2;
var x{S} >= 1 suffix stage 2;
var y{S} >= 1 suffix stage 2;
minimize o: sum{s in S} (P[s] * x[s] + s * y[s]);
