suffix stage IN;
option auxfiles rc;
set S = 3..4;
param P{S} = 1 / 2;
var x >= 1;
var y{s in S} >= s suffix stage 2;
minimize o: sum{s in S} P[s] * (x + y[s]);
