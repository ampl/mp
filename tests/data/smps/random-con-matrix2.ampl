option auxfiles rc;
suffix stage IN;
set S = 4..5;
param P{S} = 1 / card(S);
var x >= 0;
var y{S} >= 0 suffix stage 2;
minimize o: sum{s in S} P[s] * (x + 2 * y[s]);
s.t. c{s in S}: s * x + 3 * y[s] >= 6;
