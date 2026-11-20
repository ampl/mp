var x{1..5} integer >= -2 <= 2;

maximize o: -3*x[1] + 5*x[2] + 3*x[3] + 0*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: x[3] == -1;
s.t. c3: x[2] == 0;
