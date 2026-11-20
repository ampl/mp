var x{1..5} integer >= 0 <= 3;

minimize o: 2*x[1] - 2*x[2] + 4*x[3] - 5*x[4] - x[5];

s.t. c1: not (max(x[5], x[2]) <= 0);
s.t. c2:
   (if max(x[5], x[2]) <= 0 then 1 else 0)
   + (if max(x[2], x[1]) >= 1 then 1 else 0)
   + (if min(x[3], x[5]) <= 0 then 1 else 0)
   >= 3;
s.t. c3:
   (max(x[5], x[2]) <= 0)
   ==> (max(x[2], x[1]) >= 1)
   else (min(x[3], x[5]) <= 0);
s.t. c4: (min(x[3], x[5]) <= 0) <==> (max(x[5], x[2]) <= 0);
