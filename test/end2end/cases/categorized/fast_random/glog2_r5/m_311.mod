var x{1..4} integer >= 0 <= 2;

minimize o: -4*x[1] - 3*x[2] + 2*x[3] + 0*x[4];

s.t. c1: (x[4] > 0) ==> (x[1] <= 0) else (x[4] <= 3);
s.t. c2: (x[4] != -1) and (x[1] == 0) ==> (x[3] >= -1);
s.t. c3: max(x[4], x[3]) - min(x[3], x[2]) <= 2;
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) < 1;
