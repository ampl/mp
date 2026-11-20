var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] + 0*x[2] - 4*x[3];

s.t. c1: alldiff(x[3], x[1]);
s.t. c2: alldiff(x[3], x[1]);
s.t. c3: max(x[2], x[3]) - min(x[1], x[1]) < 0;
