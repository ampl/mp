var x{1..4} integer >= 0 <= 2;

maximize o: -3*x[1] - 4*x[2] + 0*x[3] - 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 2;
