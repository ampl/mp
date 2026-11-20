var x{1..4} integer >= 0 <= 3;

maximize o: -5*x[1] - 2*x[2] + 2*x[3] + 5*x[4];

s.t. c1: (if x[3] < 1 then x[1] else x[4]) != 0;
s.t. c2: not ((x[1] > 3) and (x[4] > 2));
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) <= 0;
s.t. c4: exists {i in 1..4} x[i] == 3;
