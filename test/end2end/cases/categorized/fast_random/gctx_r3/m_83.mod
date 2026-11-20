var x{1..3} integer >= -1 <= 2;

maximize o: 3*x[1] + x[2] + 6*x[3];

s.t. c1: (if x[2] < 1 then x[3] else x[3]) > 1;
s.t. c2: not ((x[2] < 1) and (min(x[3], x[2]) <= -1));
