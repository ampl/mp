var x{1..5} integer >= 0 <= 4;

maximize o: 5*x[1] - 5*x[2] - x[3] - 4*x[4] + 0*x[5];

s.t. c1: (if x[2] > 1 then x[1] else x[3]) < 4;
s.t. c2: not ((x[5] < 1) and (x[5] >= 2));
s.t. c3: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 1;
