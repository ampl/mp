var x{1..6} integer >= 0 <= 4;

minimize o: -4*x[1] - 3*x[2] + 0*x[3] + 6*x[4] + 6*x[5] + 4*x[6];

s.t. c1: not ((x[1] > 4) and (x[4] > 4));
s.t. c2: (x[5] != 1) or (x[5] == 3) or (x[5] > 3);
s.t. c3: (if x[6] != 4 then x[6] else x[3]) == 1;
s.t. c4: numberof 4 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
