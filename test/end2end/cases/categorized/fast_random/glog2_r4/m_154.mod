var x{1..5} integer >= -3 <= 1;

maximize o: 4*x[1] - 3*x[2] - 2*x[3] + 0*x[4] - 5*x[5];

s.t. c1: (x[4] != 1) or (x[4] == -1) or (x[1] < -2);
s.t. c2: exists {i in 1..5} x[i] == -3;
