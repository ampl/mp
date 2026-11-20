var x{1..4} integer >= -2 <= 3;

maximize o: x[1] + 5*x[2] - x[3] + 6*x[4];

s.t. c1: exists {i in 1..4} x[i] == -1;
