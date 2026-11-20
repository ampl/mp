var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] + 3*x[2] + x[3] - 5*x[4] + 5*x[5] + 5*x[6];

s.t. c1: max(x[6], x[3]) - min(x[1], x[4]) > 0;
s.t. c2: (x[4] != -1) or (x[3] <= 0) or (x[3] >= 0);
