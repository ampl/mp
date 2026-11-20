var x{1..6} integer >= 0 <= 1;

minimize o: -x[1] - 3*x[2] - 3*x[3] - 2*x[4] - 4*x[5] + 4*x[6];

s.t. c1: ((x[1] < 1) ==> (x[4] > 1)) and ((x[4] > 1) ==> (x[1] >= 1));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
