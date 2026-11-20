var x{1..3} integer >= -2 <= 1;

maximize o: 4*x[1] - 4*x[2] + 3*x[3];

s.t. c1: exists {i in 1..3} x[i] == 1;
