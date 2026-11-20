var x{1..5} integer >= 0 <= 4;

minimize o: 4*x[1] - 3*x[2] + 2*x[3] - 3*x[4] - x[5];

s.t. c1: exists {i in 1..5} x[i] == 0;
