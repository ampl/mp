var x{1..5} integer >= 0 <= 3;

minimize o: 3*x[1] + 0*x[2] + 0*x[3] - 2*x[4] - 4*x[5];

s.t. c1: not ((x[5] <= 1) and (x[1] >= 2));
s.t. c2: not ((x[5] <= 1) and (x[1] >= 2));
