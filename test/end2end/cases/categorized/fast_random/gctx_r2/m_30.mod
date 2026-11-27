var x{1..4} integer >= 0 <= 2;

maximize o:
   -4*x[1] - 2*x[2] + 4*x[3] - x[4] - 5*(if min(x[1], x[4]) != 2 then 1 else 0);

s.t. c1: (min(x[1], x[4]) != 2) <==> (x[1] <= 1);
s.t. c2:
   (if x[1] <= 1 then 1 else 0) + (if min(x[1], x[4]) != 2 then 1 else 0)
   + (if x[1] <= 1 then 1 else 0)
   >= 2;
s.t. c3: (x[1] <= 1) ==> ((min(x[1], x[4]) != 2) or (x[1] <= 1));
