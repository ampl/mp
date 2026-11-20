var x{1..4} integer >= -2 <= 1;

minimize o: 6*x[1] + 4*x[2] - x[3] - x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c2: (if x[4] <= -1 then x[1] else x[1]) > 1;
s.t. c3: (x[4] != 0) and (x[2] == -1) ==> (x[1] >= 0);
s.t. c4: (x[3] >= -2) and (x[2] == -1) ==> (x[3] >= -2);
