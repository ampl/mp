var x{1..4} integer >= 0 <= 3;

minimize o: 6*x[1] + 6*x[2] - 5*x[3] + 0*x[4];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c3: alldiff(x[1], x[3]);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) != 3;
