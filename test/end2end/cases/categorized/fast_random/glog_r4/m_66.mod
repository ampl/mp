var x{1..4} integer >= -2 <= 1;

maximize o: 5*x[1] + 4*x[2] + 5*x[3] + x[4];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4]) >= 1;
