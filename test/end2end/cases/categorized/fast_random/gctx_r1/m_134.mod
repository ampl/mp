var x{1..4} integer >= -1 <= 1;

maximize o:
   -4*x[1] + 6*x[2] + 2*x[3] + 4*x[4] + 5*(if x[3] <= -1 then 1 else 0);

s.t. c1:
   (if x[3] <= -1 then 1 else 0) + (if x[3] <= -1 then 1 else 0)
   + (if x[3] <= -1 then 1 else 0)
   >= 1;
s.t. c2: not ((x[3] <= -1) and (min(x[4], x[3]) >= 1));
s.t. c3: (if x[3] <= -1 then x[2] else x[4]) > 0;
