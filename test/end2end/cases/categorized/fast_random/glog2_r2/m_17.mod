var x{1..4} integer >= 0 <= 2;

maximize o: 1000*x[1] + 51*x[2] + 51*x[3] - 333*x[4];

s.t. c1: max(x[3], x[1]) - min(x[4], x[3]) < 1;
s.t. c2: (if x[3] <= 1 then x[1] else x[3]) <= 1;
s.t. c3: (if x[4] <= 2 then x[2] else x[2]) < 0;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4]) < 0;
