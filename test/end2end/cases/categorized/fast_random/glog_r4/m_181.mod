var x{1..3} integer >= -2 <= 2;

minimize o: 6*x[1] + 5*x[2] - 5*x[3];

s.t. c1: exists {i in 1..3} x[i] == 1;
