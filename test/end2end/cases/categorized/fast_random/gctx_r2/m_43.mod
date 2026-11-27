var x{1..5} integer >= -1 <= 2;

minimize o:
   -4*x[1] + 2*x[2] + x[3] - 4*x[4] + 3*x[5]
   - 4*(if min(x[4], x[4]) != 0 then 1 else 0);

s.t. c1: (x[5] > 1) ==> (min(x[4], x[4]) != 0);
s.t. c2:
   (if min(x[4], x[4]) != 0 then 1 else 0)
   + (if min(x[4], x[4]) != 0 then 1 else 0)
   + (if x[2] + x[2] <= -1 then 1 else 0)
   >= 3;
s.t. c3: (min(x[4], x[4]) != 0) <==> (x[1] == x[4]);
