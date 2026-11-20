var x{1..4} integer >= 0 <= 3;

minimize o: x[1] - 3*x[2] + 4*x[3] - 2*x[4] + (if x[4] != x[2] then 1 else 0);

s.t. c1: (x[4] != x[2]) <==> (x[4] != x[2]);
s.t. c2: (x[4] != x[2]) or (x[4] != x[2]);
s.t. c3: (x[4] != x[2]) ==> (x[4] != x[2]) else (x[4] != x[2]);
