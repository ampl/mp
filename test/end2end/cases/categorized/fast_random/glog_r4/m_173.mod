var x{1..6} integer >= -2 <= 1;

minimize o: x[1] + 6*x[2] - 5*x[3] + 3*x[4] - 2*x[5] - 4*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[6] > 1) or (x[3] >= -2) or (x[5] != 0);
s.t. c3: (x[3] != 0) ==> (x[2] != -1);
