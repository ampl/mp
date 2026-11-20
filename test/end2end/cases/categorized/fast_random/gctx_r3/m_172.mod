var x{1..4} integer >= -1 <= 1;

minimize o:
   2*x[1] - 4*x[2] - 4*x[3] - x[4] + 5*(if abs(x[4] - x[4]) > 2 then 1 else 0);

s.t. c1: (x[2] < 1) ==> (x[2] < 1) else (x[2] < 1);
s.t. c2:
   (if x[2] < 1 then 1 else 0) + (if min(x[4], x[4]) >= 0 then 1 else 0)
   + (if abs(x[4] - x[4]) > 2 then 1 else 0)
   >= 3;
