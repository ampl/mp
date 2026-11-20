var x{1..6} integer >= 0 <= 4;

maximize o: x[1] + 5*x[2] + 5*x[3] + 6*x[4] + 3*x[5] + 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 3;
