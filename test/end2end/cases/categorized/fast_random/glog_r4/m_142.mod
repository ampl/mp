var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] + 2*x[2] + 0*x[3] + 4*x[4];

s.t. c1: (x[4] >= 4) ==> (x[4] <= 1);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) == 2;
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4]) < 2;
