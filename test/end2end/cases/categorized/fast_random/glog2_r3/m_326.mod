var x{1..4} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: not ((x[2] >= 0) and (x[3] <= 2));
s.t. c2: (x[1] != 0) or (x[3] > 2) or (x[4] <= 0);
