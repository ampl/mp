var x{1..4} integer >= -2 <= 0;

maximize o: 4*x[1] - x[2] + x[3] + 6*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) <= 3;
s.t. c2: not ((x[2] <= -1) and (x[3] <= 0));
s.t. c3: not ((x[1] <= -1) and (x[1] < 0));
s.t. c4: not ((x[3] <= -1) and (x[1] == -1));
s.t. c5: x[1] == -1;
