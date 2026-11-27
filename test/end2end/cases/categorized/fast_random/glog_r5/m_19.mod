var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + x[2] - 2*x[3] - 2*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) <= 0;
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: alldiff(x[2], x[1], x[4], x[3]);
