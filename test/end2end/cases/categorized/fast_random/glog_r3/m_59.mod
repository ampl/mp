var x{1..4} integer >= 0 <= 3;

minimize o: 5*x[1] + x[2] - 3*x[3] + 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
