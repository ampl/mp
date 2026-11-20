var x{1..4} integer >= -1 <= 2;

minimize o: 5*x[1] + 4*x[2] + 3*x[3] + 0*x[4];

s.t. c1: (x[2] >= 0) <==> (min(x[4], x[1]) > 2);
s.t. c2:
   (if min(x[4], x[1]) > 2 then x[3] else 0) + (if x[1] != 1 then x[2] else 0)
   >= 1;
s.t. c3: (min(x[4], x[1]) > 2) <==> (x[2] >= 0);
