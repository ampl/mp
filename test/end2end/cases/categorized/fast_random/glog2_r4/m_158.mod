var x{1..4} integer >= 0 <= 1;

minimize o: -4*x[1] + 3*x[2] + 2*x[3] - 2*x[4];

s.t. c1: (x[4] < 1) <==> (x[2] + x[3] <= 0);
s.t. c2: max(x[3], x[4]) - min(x[2], x[4]) <= 0;
