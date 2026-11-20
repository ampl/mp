var x{1..4} integer >= 0 <= 4;

minimize o: 2*x[1] + 3*x[2] - x[3] - 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
