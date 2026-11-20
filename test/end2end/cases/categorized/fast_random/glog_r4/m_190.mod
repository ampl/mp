var x{1..6} integer >= 0 <= 4;

minimize o: x[1] + 5*x[2] + 3*x[3] + 3*x[4] + x[5] - 5*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 3;
s.t. c2: (x[4] <= 4) or (x[5] == 4) or (x[3] == 3);
s.t. c3: (if x[2] <= 0 then x[4] else x[4]) >= 3;
