var x{1..5} integer >= 0 <= 4;

minimize o: 2*x[1] - x[2] - 4*x[3] - 3*x[4] + 2*x[5];

s.t. c1: not ((x[5] >= 1) and (x[5] <= 0));
s.t. c2: (x[4] != 1) <==> (x[4] + x[2] > 3);
