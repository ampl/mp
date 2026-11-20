var x{1..4} integer >= -2 <= 1;

maximize o: -4*x[1] - x[2] + 2*x[3] - x[4];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4]) != 0;
s.t. c2: exists {i in 1..4} x[i] == -2;
s.t. c3: numberof -2 in (x[1], x[2], x[3], x[4]) != 0;
s.t. c4: exists {i in 1..4} x[i] == -2;
