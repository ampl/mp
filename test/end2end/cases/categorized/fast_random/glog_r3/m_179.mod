var x{1..5} integer >= 0 <= 4;

maximize o: 3*x[1] + 2*x[2] + 5*x[3] + 3*x[4] + 2*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
