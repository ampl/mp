var x{1..4} integer >= -2 <= 0;

maximize o: -2*x[1] - 3*x[2] - 4*x[3] - 5*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) >= 3;
