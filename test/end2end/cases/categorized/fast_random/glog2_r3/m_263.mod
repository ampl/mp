var x{1..6} integer >= 0 <= 4;

maximize o: 4*x[1] - 2*x[2] + 4*x[3] + 0*x[4] + 2*x[5] - 5*x[6];

s.t. c1: not ((x[6] > 4) and (x[2] != 4));
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 0;
s.t. c3: x[2] == 4;
