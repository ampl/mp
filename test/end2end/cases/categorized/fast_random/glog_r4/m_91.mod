var x{1..4} integer >= 0 <= 2;

maximize o: 2*x[1] + 0*x[2] - x[3] + 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
