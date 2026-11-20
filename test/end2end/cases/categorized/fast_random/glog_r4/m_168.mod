var x{1..6} integer >= -2 <= 1;

maximize o: -4*x[1] + 0*x[2] - x[3] - 4*x[4] + 2*x[5] + 3*x[6];

s.t. c1: (x[1] >= -1) ==> (x[2] != 0);
s.t. c2: abs(x[3] - x[5]) >= 1;
s.t. c3: exists {i in 1..6} x[i] == 1;
