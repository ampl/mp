var x{1..4} integer >= 0 <= 2;

maximize o:
   -5*x[1] + x[2] + 6*x[3] - 5*x[4] + 2*(if min(x[1], x[2]) < 0 then 1 else 0);

s.t. c1: (if min(x[1], x[2]) < 0 then x[1] else x[1]) >= 0;
s.t. c2:
   (if min(x[1], x[2]) < 0 then 1 else 0)
   + (if min(x[1], x[2]) < 0 then 1 else 0)
   + (if min(x[1], x[2]) < 0 then 1 else 0)
   >= 2;
s.t. c3: (if min(x[1], x[2]) < 0 then x[1] else x[1]) > 2;
