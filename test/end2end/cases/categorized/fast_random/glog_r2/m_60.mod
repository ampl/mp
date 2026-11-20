var x{1..4} integer >= -2 <= 2;

maximize o: -3*x[1] + 6*x[2] + 6*x[3] + 6*x[4];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) == 0;
s.t. c2: alldiff(x[3], x[4]);
