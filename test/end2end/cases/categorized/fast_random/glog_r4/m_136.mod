var x{1..5} integer >= 0 <= 2;

minimize o: -4*x[1] - 4*x[2] + 2*x[3] + 4*x[4] + x[5];

s.t. c1: (x[3] <= 0) <==> (x[5] + x[2] > 0);
s.t. c2: abs(x[4] - x[4]) > 1;
s.t. c3: max(x[4], x[1]) - min(x[1], x[2]) != 2;
