option auxfiles rc;
suffix stage IN;
set S = 0..1;
param P{S} = 1 / 2;
var x >= 1;
var z >= 2;
var y{S} >= 3 suffix stage 2;
minimize o: sum{s in S} P[s] * (x + z + y[s]);
# c[0] should be eliminated by presolve
s.t. c{s in S}: s * x - s * y[s] >= 0;
