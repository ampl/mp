var x{1..3} integer >= 0 <= 3;

minimize o: 2000*x[1] - 1000*x[2] - 999*x[3];

s.t. c1: not ((x[1] <= 1) and (x[2] <= 0));
s.t. c2: max(x[3], x[2]) - min(x[3], x[3]) != 3;
