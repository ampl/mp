var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (if x[3] > 0 then x[1] else x[4]) <= 0;
s.t. c2: (if x[4] != 1 then x[4] else x[2]) == 2;
