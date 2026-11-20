var x{1..6} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: not ((x[5] >= 0) and (x[6] == 2));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (x[4] != 1) or (x[3] >= 0) or (x[2] <= 2);
