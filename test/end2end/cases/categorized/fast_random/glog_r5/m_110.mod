var x{1..5} integer >= 0 <= 4;

minimize o: -5*x[1] - 5*x[2] + 0*x[3] + 2*x[4] + 2*x[5];

s.t. c1: exists {i in 1..5} x[i] == 4;
s.t. c2: not ((x[4] <= 4) and (x[4] != 1));
