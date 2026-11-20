var x{1..4} integer >= -2 <= 2;

minimize o: -999*x[1] + 17*x[2] + 2000*x[3] - 2000*x[4];

s.t. c1: abs(x[3] - x[4]) >= 4;
s.t. c2: alldiff(x[3], x[4], x[1], x[2]);
s.t. c3: abs(x[2] - x[4]) >= 2;
