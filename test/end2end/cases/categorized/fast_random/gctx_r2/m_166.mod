var x{1..4} integer >= 0 <= 3;

minimize o: 4*x[1] - 2*x[2] + 0*x[3] + 4*x[4];

s.t. c1: (if x[2] > 1 then x[2] else 0) + (if x[2] > 1 then x[1] else 0) > 2;
s.t. c2: (x[3] == x[3]) or (x[3] == x[3]);
s.t. c3: not (x[3] == x[3]);
