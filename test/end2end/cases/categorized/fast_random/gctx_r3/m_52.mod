var x{1..5} integer >= -1 <= 1;

maximize o: 4*x[1] + x[2] + 5*x[3] + 3*x[4] + 0*x[5];

s.t. c1: (min(x[4], x[1]) > 0) ==> (x[2] <= 0) else (min(x[4], x[1]) > 0);
s.t. c2: not ((min(x[4], x[1]) > 0) and (min(x[4], x[1]) > 0));
s.t. c3:
   (if min(x[4], x[1]) > 0 then x[5] else 0)
   + (if abs(x[2] - x[3]) < 1 then x[2] else 0)
   <= 0;
s.t. c4: (min(x[4], x[1]) > 0) or (min(x[4], x[1]) < 0);
