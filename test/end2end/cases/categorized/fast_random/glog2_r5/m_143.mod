var x{1..4} integer >= 0 <= 4;

maximize o: -4*x[1] - 2*x[2] + 2*x[3] - 4*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) != 0;
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: x[3] == 0;
