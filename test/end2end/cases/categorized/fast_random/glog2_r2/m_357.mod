var x{1..6} integer >= 0 <= 4;

maximize o: 17*x[1] + 3000*x[2] + 51*x[3] + 2000*x[4] + 51*x[5] - 1000*x[6];

s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
s.t. c2: abs(x[4] - x[4]) == 4;
