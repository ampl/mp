var x{1..4} integer >= -2 <= 3;

minimize o: -4*x[1] + 6*x[2] - 3*x[3] + 2*x[4];

s.t. c1: max(x[1], x[3]) - min(x[4], x[1]) <= 1;
s.t. c2: alldiff(x[3], x[2]);
s.t. c3: abs(x[4] - x[2]) < 4;
