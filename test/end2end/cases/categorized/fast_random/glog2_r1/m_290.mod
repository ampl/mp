var x{1..5} integer >= 0 <= 2;

maximize o: -4*x[1] - 5*x[2] - 5*x[3] + x[4] - x[5];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
s.t. c2: alldiff(x[2], x[4], x[1], x[5]);
s.t. c3: max(x[4], x[4]) - min(x[3], x[2]) != 2;
