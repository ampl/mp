var x{1..4} integer >= 0 <= 4;

minimize o: 5*x[1] + 2*x[2] - 2*x[3] + 4*x[4];

s.t. c1: numberof 4 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c2: alldiff(x[3], x[1], x[2], x[4]);
s.t. c3: numberof 4 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c4: alldiff(x[3], x[1], x[2], x[4]);
