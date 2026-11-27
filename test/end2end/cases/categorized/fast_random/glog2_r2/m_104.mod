var x{1..5} integer >= -3 <= 1;

maximize o: -2*x[1] - 5*x[2] + 0*x[3] - x[4] + 6*x[5];

s.t. c1: exists {i in 1..5} x[i] == -3;
s.t. c2: not ((x[4] >= 0) and (x[3] == 0));
