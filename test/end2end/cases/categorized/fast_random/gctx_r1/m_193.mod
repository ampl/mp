var x{1..4} integer >= -1 <= 1;

maximize o:
   -3*x[1] + 3*x[2] + 0*x[3] + 5*x[4]
   + 5*(if min(x[3], x[4]) != -1 then 1 else 0);

s.t. c1: (min(x[3], x[4]) != -1) ==> (min(x[3], x[4]) != -1);
s.t. c2:
   (if abs(x[1] - x[1]) >= 0 then x[4] else 0)
   + (if abs(x[1] - x[1]) >= 0 then x[1] else 0)
   < 1;
