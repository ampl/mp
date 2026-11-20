var x{1..5} integer >= -1 <= 2;

minimize o:
   4*x[1] - 3*x[2] + 2*x[3] + x[4] - 2*x[5]
   + 3*(if max(x[3], x[5]) <= 2 then 1 else 0);

s.t. c1:
   (if x[4] >= x[4] then x[1] else 0)
   + (if min(x[3], x[4]) > 0 then x[2] else 0)
   > -1;
s.t. c2: (x[4] >= x[4]) <==> (min(x[3], x[4]) > 0);
