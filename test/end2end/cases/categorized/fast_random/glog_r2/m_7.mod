var x{1..4} integer >= -2 <= 1;

maximize o: x[1] + 5*x[2] - 5*x[3] - 4*x[4];

s.t. c1: (x[4] >= 1) or (x[4] >= -2) or (x[1] <= -1);
s.t. c2: abs(x[2] - x[4]) == 2;
s.t. c3: exists {i in 1..4} x[i] == -1;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) != 1;
