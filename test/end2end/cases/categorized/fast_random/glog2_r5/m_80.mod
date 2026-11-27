var x{1..6} integer >= 0 <= 2;

minimize o: 3*x[1] + 6*x[2] - 4*x[3] + x[4] + x[5] + 0*x[6];

s.t. c1: (x[5] > 0) <==> (x[5] + x[2] > 2);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
