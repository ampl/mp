suffix stage IN;
option auxfiles rc;
set S = {3};
param P{S} = 1;
var x >= 1;
var y{s in S} >= 2 suffix stage 2;
minimize o: sum{s in S} P[s] * (x + y[s]);
