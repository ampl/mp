var x{1..5} integer >= 0 <= 2;

minimize o: 3*x[1] - 3*x[2] - x[3] - x[4] - 4*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: max(x[2], x[5]) - min(x[5], x[4]) != 2;
s.t. c3: not ((x[2] < 2) and (x[4] != 0));
