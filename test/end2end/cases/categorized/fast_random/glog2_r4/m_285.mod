var x{1..4} integer >= 0 <= 2;

minimize o: -3*x[1] + 2*x[2] - 4*x[3] - 2*x[4];

s.t. c1: abs(x[3] - x[4]) >= 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c3: exists {i in 1..4} x[i] == 2;
s.t. c4: exists {i in 1..4} x[i] == 1;
s.t. c5: not ((x[2] >= 0) and (x[4] < 0));
s.t. c6: (if x[4] > 0 then x[1] else x[4]) > 2;
