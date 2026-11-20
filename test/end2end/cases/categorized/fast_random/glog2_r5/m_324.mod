var x{1..3} integer >= -3 <= 1;

minimize o: 0*x[1] + 3*x[2] + 0*x[3];

s.t. c1: alldiff(x[3], x[2], x[1]);
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: not ((x[2] != 1) and (x[3] >= 0));
s.t. c4: numberof -3 in (x[1], x[2], x[3]) != 0;
