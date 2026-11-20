var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] + 4*x[2] + 0*x[3] - x[4];

s.t. c1: not (x[2] == x[3]);
s.t. c2:
   (if x[2] == x[3] then x[4] else 0) + (if x[2] == x[3] then x[2] else 0) > 2;
s.t. c3: (x[2] == x[3]) ==> (min(x[2], x[1]) <= 1);
