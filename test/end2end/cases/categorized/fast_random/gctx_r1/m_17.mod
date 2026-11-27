var x{1..5} integer >= 0 <= 3;

maximize o: 6*x[1] - 2*x[2] + 2*x[3] - 3*x[4] - x[5];

s.t. c1: not ((x[2] >= 1) and (x[2] != x[4]));
s.t. c2: (if min(x[1], x[4]) <= 1 then x[5] else x[1]) >= 2;
s.t. c3: (min(x[1], x[4]) <= 1) and ((min(x[1], x[4]) <= 1) or not (x[2] >= 1));
