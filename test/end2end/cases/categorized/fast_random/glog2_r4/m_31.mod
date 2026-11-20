var x{1..5} integer >= 0 <= 4;

minimize o: 5*x[1] + x[2] - 4*x[3] + x[4] - 5*x[5];

s.t. c1: not ((x[2] >= -1) and (x[4] >= 4));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
