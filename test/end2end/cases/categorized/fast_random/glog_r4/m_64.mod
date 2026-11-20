var x{1..6} integer >= -2 <= 0;

minimize o: -4*x[1] + 6*x[2] + 4*x[3] - 2*x[4] - 2*x[5] - 3*x[6];

s.t. c1: alldiff(x[2], x[4], x[5]);
s.t. c2: abs(x[1] - x[1]) == 2;
