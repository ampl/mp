var x{1..4} integer >= -1 <= 1;

minimize o: -4*x[1] + x[2] + 3*x[3] + 6*x[4];

s.t. c1: not (x[3] + x[3] == -1);
s.t. c2: (x[4] <= x[3]) <==> (abs(x[1] - x[1]) > 2);
