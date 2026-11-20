var x{1..3} integer >= 0 <= 3;

maximize o: -3*x[1] + 6*x[2] - 5*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) >= 1;
