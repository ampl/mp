var x{1..4} integer >= -2 <= 0;

minimize o: -4*x[1] + x[2] + 4*x[3] + 0*x[4];

s.t. c1: max(x[1], x[3]) - min(x[4], x[3]) >= -1;
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c4: not ((x[1] > 0) and (x[1] >= -1));
s.t. c5: x[2] == -1;
s.t. c6: x[3] == -2;
