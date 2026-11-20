var x{1..4} integer >= 0 <= 2;

minimize o: 0*x[1] - 3*x[2] + 5*x[3] + x[4] + 3*(if x[1] > 0 then 1 else 0);

s.t. c1: not (x[1] != x[1]);
s.t. c2: (x[1] != x[1]) <==> (x[1] > 0);
s.t. c3: not ((x[1] > 0) and (x[1] != x[1]));
