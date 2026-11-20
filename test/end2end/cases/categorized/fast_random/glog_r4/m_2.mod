var x{1..6} integer >= -2 <= 0;

minimize o: 2*x[1] - 4*x[2] - 4*x[3] + 4*x[4] + x[5] + 0*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
