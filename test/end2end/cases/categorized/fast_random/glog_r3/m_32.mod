var x{1..4} integer >= 0 <= 3;

maximize o: 5*x[1] - 4*x[2] - 4*x[3] + x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
