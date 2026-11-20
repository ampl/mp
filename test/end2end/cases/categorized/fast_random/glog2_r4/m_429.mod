var x{1..6} integer >= 0 <= 4;

maximize o: -5*x[1] + 2*x[2] - 4*x[3] + 4*x[4] - 2*x[5] - 3*x[6];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 1;
s.t. c2: (if x[4] > 4 then x[2] else x[1]) < 3;
s.t. c3: abs(x[3] - x[4]) != 1;
s.t. c4: not ((x[6] == 3) and (x[2] == 3));
s.t. c5: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 1;
s.t. c6: not ((x[6] == 3) and (x[2] == 3));
