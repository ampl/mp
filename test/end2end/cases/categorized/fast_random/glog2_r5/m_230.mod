var x{1..4} integer >= 0 <= 3;

maximize o: 5*x[1] + 0*x[2] + x[3] - 5*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
