var x{1..4} integer >= 0 <= 2;

minimize o: x[1] + 6*x[2] + 2*x[3] + 4*x[4];

s.t. c1: not (x[4] >= 2);
s.t. c2: (x[1] + x[3] < 2) <==> (x[4] >= 2);
s.t. c3: not (x[4] <= 0);
