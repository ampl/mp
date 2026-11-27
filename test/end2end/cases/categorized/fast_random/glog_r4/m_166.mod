var x{1..4} integer >= -2 <= 1;

maximize o: -5*x[1] + 4*x[2] + 0*x[3] + 6*x[4];

s.t. c1: (if x[1] >= 1 then x[1] else x[1]) == 1;
s.t. c2: max(x[4], x[4]) - min(x[2], x[1]) == -1;
s.t. c3: (if x[2] >= -2 then x[1] else x[4]) >= -1;
