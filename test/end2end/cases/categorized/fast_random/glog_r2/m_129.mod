var x{1..6} integer >= 0 <= 2;

minimize o: 0*x[1] - 5*x[2] + 2*x[3] + 3*x[4] + 5*x[5] + 5*x[6];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 0;
s.t. c2: abs(x[5] - x[2]) <= 1;
s.t. c3: alldiff(x[6], x[4]);
