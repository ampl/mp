var x{1..6} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] - 4*x[3] + x[4] - 2*x[5] - 5*x[6];

s.t. c1: not ((x[3] >= 0) and (x[6] >= 0));
s.t. c2: not ((x[3] <= 0) and (x[2] >= 1));
s.t. c3: exists {i in 1..6} x[i] == 3;
s.t. c4: not ((x[3] >= 0) and (x[6] >= 0));
s.t. c5: exists {i in 1..6} x[i] == 3;
