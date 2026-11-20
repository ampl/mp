var x{1..6} integer >= -2 <= 1;

minimize o: x[1] - 2*x[2] - 4*x[3] + x[4] + 3*x[5] + 4*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[3]);
