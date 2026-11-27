var x{1..6} integer >= 0 <= 4;

minimize o: 5*x[1] - 5*x[2] + 0*x[3] + 5*x[4] - 3*x[5] - 3*x[6];

s.t. c1: abs(x[1] - x[3]) == 2;
s.t. c2: not ((x[1] > 2) and (x[4] > 3));
s.t. c3: exists {i in 1..6} x[i] == 4;
