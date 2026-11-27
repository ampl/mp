var x{1..5} integer >= -1 <= 2;

minimize o:
   -x[1] - 5*x[2] + x[3] - 4*x[4] - 4*x[5]
   + 6*(if min(x[1], x[1]) != 0 then 1 else 0);

s.t. c1:
   (if min(x[5], x[1]) != 0 then x[1] else 0)
   + (if min(x[1], x[1]) != 0 then x[4] else 0)
   < 1;
s.t. c2: (min(x[1], x[1]) != 0) <==> (min(x[5], x[1]) != 0);
s.t. c3:
   (min(x[5], x[1]) != 0)
   and ((min(x[5], x[1]) != 0) or not (min(x[5], x[1]) != 0));
s.t. c4:
   (min(x[1], x[1]) != 0)
   ==> ((min(x[5], x[1]) != 0) or (min(x[5], x[1]) != 0));
