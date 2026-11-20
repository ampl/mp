var x{1..6} integer >= 0 <= 2;

minimize o: -x[1] - 4*x[2] + x[3] - x[4] + 0*x[5] + 4*x[6];

s.t. c1: (x[3] < 2) <==> (x[2] + x[6] <= 0);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 0;
