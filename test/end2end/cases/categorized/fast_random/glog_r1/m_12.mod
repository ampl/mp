var x{1..5} integer >= 0 <= 3;

maximize o: 0*x[1] + x[2] + x[3] + 5*x[4] + 3*x[5];

s.t. c1: not ((x[3] <= 0) and (x[3] == 1));
s.t. c2: exists {i in 1..5} x[i] == 2;
