var x{1..5} integer >= 0 <= 4;

maximize o: 5*x[1] + 0*x[2] + 4*x[3] + 3*x[4] - 4*x[5];

s.t. c1: max(x[4], x[4]) - min(x[4], x[4]) > 4;
s.t. c2: max(x[5], x[5]) - min(x[2], x[5]) >= 2;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
