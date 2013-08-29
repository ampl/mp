option auxfiles rc;
suffix stage IN;
set S = 0..1;
param P{S} = 1 / 2;
var x >= 0;
var y >= 0;
var z{S} >= 0 suffix stage 2;
minimize o: sum{s in S} P[s] * (x + y + z[s]);
s.t. c{s in S}: x + y + s * z[s] >= 5;
