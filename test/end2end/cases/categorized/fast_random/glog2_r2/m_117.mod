var x{1..4} integer >= -2 <= 1;

minimize o: -2000*x[1] - 999*x[2] + 34*x[3] - 999*x[4];

s.t. c1: (x[2] < 1) and (x[2] > -2) ==> (x[3] >= 1);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) > 2;
