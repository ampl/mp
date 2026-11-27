var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] - 4*x[2] + 5*x[3] - 4*x[4] + 2*x[5] + 2*x[6];

s.t. c1: (if x[5] >= 1 then x[2] else x[2]) == 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 1;
s.t. c3: (if x[5] >= 1 then x[2] else x[2]) == 1;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 1;
