var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] + 5*x[2] - x[3] + 6*x[4] + 6*x[5] - 4*x[6];

s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
s.t. c2: max(x[4], x[6]) - min(x[3], x[2]) > 1;
s.t. c3: x[6] == 4;
s.t. c4: x[1] == 2;
