var x{1..4} integer >= -1 <= 2;

minimize o: -4*x[1] - 5*x[2] - x[3] + 3*x[4];

s.t. c1: (x[3] > x[1]) <==> (x[3] > x[1]);
s.t. c2: not (x[3] > x[1]);
