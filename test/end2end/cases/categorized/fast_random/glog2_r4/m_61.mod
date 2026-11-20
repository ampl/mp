var x{1..3} integer >= 0 <= 3;

maximize o: 3*x[1] - 3*x[2] - 5*x[3];

s.t. c1: not ((x[2] < 0) and (x[1] == -1));
s.t. c2: (if x[2] != 0 then x[1] else x[3]) >= 4;
s.t. c3: max(x[2], x[1]) - min(x[3], x[3]) >= 3;
