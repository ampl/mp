var x{1..6} integer >= 0 <= 2;

minimize o: 3*x[1] - 3*x[2] + 4*x[3] - 5*x[4] - 4*x[5] + 5*x[6];

s.t. c1: (if x[1] > 0 then x[6] else x[2]) > 0;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
