var x{1..5} integer >= -2 <= 0;

minimize o: -2*x[1] + 0*x[2] + 6*x[3] - x[4] - 4*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
