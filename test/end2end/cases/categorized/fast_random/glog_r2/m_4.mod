var x{1..6} integer >= 0 <= 2;

minimize o: 4*x[1] + 0*x[2] + 2*x[3] + 3*x[4] - 2*x[5] + 5*x[6];

s.t. c1: (if x[5] == 0 then x[3] else x[4]) == 0;
s.t. c2: max(x[5], x[2]) - min(x[5], x[2]) >= 2;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
