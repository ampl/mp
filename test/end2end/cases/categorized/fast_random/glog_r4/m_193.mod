var x{1..3} integer >= 0 <= 4;

maximize o: 3*x[1] - 3*x[2] - x[3];

s.t. c1: (if x[2] >= 1 then x[3] else x[2]) < 4;
s.t. c2: numberof 2 in (x[1], x[2], x[3]) == 0;
