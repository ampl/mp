var x{1..6} integer >= 0 <= 3;

maximize o: 4*x[1] + 3*x[2] + 6*x[3] + 0*x[4] - 2*x[5] + 0*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
