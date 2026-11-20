var x{1..4} integer >= 0 <= 2;

maximize o:
   -2*x[1] + 3*x[2] - 5*x[3] + x[4] + 3*(if min(x[2], x[1]) <= 0 then 1 else 0);

s.t. c1:
   (if min(x[3], x[2]) != 2 then 1 else 0)
   + (if min(x[3], x[2]) != 2 then 1 else 0)
   + (if min(x[2], x[1]) <= 0 then 1 else 0)
   >= 3;
s.t. c2:
   (if min(x[3], x[2]) != 2 then x[2] else 0)
   + (if min(x[2], x[1]) <= 0 then x[3] else 0)
   >= 1;
s.t. c3: (min(x[2], x[1]) <= 0) <==> (min(x[2], x[1]) <= 0);
s.t. c4:
   (min(x[2], x[1]) <= 0) ==> (min(x[2], x[1]) < 2) else (min(x[2], x[1]) < 2);
