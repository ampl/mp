var x{1..7} integer >= 0 <= 1;

maximize o: -5*x[1] + 4*x[2] + 0*x[3] - 3*x[4] + 5*x[5] + x[6] - 2*x[7];

s.t. c1: forall {i in 1..6} x[i] != x[i+1];
s.t. c2: (x[7] <= 1) or (x[6] <= 1) or (x[4] < 1);
