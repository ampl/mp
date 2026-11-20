var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] + 4*x[3] + 6*x[4];

s.t. c1: (x[2] <= 2) <==> (x[4] + x[4] < 3);
s.t. c2: alldiff(x[1], x[4], x[3], x[2]);
