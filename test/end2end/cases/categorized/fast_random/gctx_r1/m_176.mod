var x{1..4} integer >= -1 <= 2;

minimize o:
   x[1] - 4*x[2] + 5*x[3] + 4*x[4] - 5*(if max(x[1], x[2]) == 2 then 1 else 0);

s.t. c1:
   (if min(x[1], x[3]) != 0 then 1 else 0)
   + (if min(x[1], x[3]) != 0 then 1 else 0)
   + (if min(x[1], x[3]) != 0 then 1 else 0)
   >= 2;
s.t. c2: (max(x[1], x[2]) == 2) ==> (max(x[2], x[4]) > 2);
s.t. c3:
   (min(x[1], x[3]) != 0) ==> (max(x[1], x[2]) == 2) else (max(x[2], x[4]) > 2);
s.t. c4:
   (min(x[1], x[3]) != 0) ==> (max(x[2], x[4]) > 2) else (max(x[1], x[2]) == 2);
