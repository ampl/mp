var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] + 5*x[2] - 3*x[3] + x[4];

s.t. c1: (min(x[3], x[4]) > 2) <==> (x[3] != x[3]);
s.t. c2: not ((x[3] != x[3]) and (x[3] != x[3]));
s.t. c3: (x[3] != x[3]) <==> (x[3] != x[3]);
