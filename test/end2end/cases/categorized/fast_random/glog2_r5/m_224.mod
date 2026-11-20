var x{1..5} integer >= -3 <= 2;

minimize o: -4*x[1] - 4*x[2] + 0*x[3] + 4*x[4] - x[5];

s.t. c1: (x[4] <= -1) <==> (x[5] + x[3] <= 2);
s.t. c2: alldiff(x[4], x[1], x[5]);
