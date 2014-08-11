suffix stage IN;
option auxfiles rc;
set S = 3..4;
param P{S} = 1 / 2;
var y{s in S} >= 2 suffix stage 2;
var x >= 1;
minimize o: sum{s in S} P[s] * (x + y[s]);
s.t. c{s in S}: 1.5 <= x - y[s];
