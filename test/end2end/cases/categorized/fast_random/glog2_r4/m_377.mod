var x{1..3} integer >= 0 <= 4;

maximize o: 34*x[1] + 1000*x[2] + 51*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) > 3;
s.t. c2: (if x[2] < 2 then x[2] else x[1]) != 4;
s.t. c3: ((if x[3] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c4: numberof 1 in (x[1], x[2], x[3]) == 1;
