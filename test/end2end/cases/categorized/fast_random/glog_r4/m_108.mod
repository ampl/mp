var x{1..6} integer >= -2 <= 1;

maximize o: 0*x[1] - 5*x[2] - 2*x[3] + 3*x[4] + 4*x[5] + 0*x[6];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 2;
