var x{1..6} integer >= 0 <= 3;

maximize o: -3*x[1] - x[2] - 4*x[3] + 0*x[4] - 3*x[5] + 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 3;
s.t. c2: abs(x[5] - x[6]) == 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
