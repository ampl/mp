var x{1..3} integer >= 0 <= 3;

minimize o: 4*x[1] + 2*x[2] - 5*x[3];

s.t. c1: alldiff(x[2], x[3]);
s.t. c2: abs(x[3] - x[3]) > 2;
