var x{1..4} integer >= 0 <= 3;

minimize o:
   -2*x[1] + x[2] - 2*x[3] + x[4] + 2*(if min(x[2], x[3]) < 2 then 1 else 0);

s.t. c1: (max(x[1], x[2]) < 1) ==> ((x[1] != 2) or (max(x[1], x[2]) < 1));
s.t. c2:
   (if min(x[2], x[3]) < 2 then 1 else 0) + (if x[1] != 2 then 1 else 0)
   + (if min(x[2], x[3]) < 2 then 1 else 0)
   >= 2;
s.t. c3: (min(x[2], x[3]) < 2) <==> (min(x[2], x[3]) < 2);
