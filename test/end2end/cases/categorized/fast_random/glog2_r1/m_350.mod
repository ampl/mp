var x{1..6} integer >= 0 <= 3;

minimize o: -4*x[1] - 3*x[2] - 5*x[3] + 4*x[4] - x[5] + 4*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
