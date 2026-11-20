var x{1..3} integer >= 0 <= 2;

minimize o: 5*x[1] + x[2] - 3*x[3] + (if min(x[2], x[1]) <= 2 then 1 else 0);

s.t. c1: (min(x[2], x[1]) <= 2) and ((min(x[2], x[1]) <= 2) or not (x[3] != 1));
s.t. c2: (x[3] != 1) <==> (x[3] != 1);
s.t. c3: not (x[3] != 1);
