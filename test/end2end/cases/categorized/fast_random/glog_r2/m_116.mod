var x{1..4} integer >= 0 <= 2;

minimize o: 3*x[1] - 4*x[2] + 0*x[3] - 2*x[4];

s.t. c1: (x[4] < 2) ==> (x[4] == 2) else (x[2] < 2);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) >= 2;
