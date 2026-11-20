var x{1..5} integer >= 0 <= 3;

maximize o:
   5*x[1] + x[2] + 3*x[3] - 3*x[4] + 0*x[5]
   + 4*(if max(x[2], x[3]) == 1 then 1 else 0);

s.t. c1:
   (if max(x[2], x[3]) == 1 then 1 else 0)
   + (if min(x[3], x[4]) > 0 then 1 else 0) + (if x[5] <= x[4] then 1 else 0)
   >= 2;
s.t. c2: not ((min(x[5], x[3]) > 1) and (x[5] <= x[4]));
s.t. c3:
   (max(x[2], x[3]) == 1) and ((x[5] <= x[4]) or not (min(x[5], x[3]) > 1));
s.t. c4:
   (if min(x[3], x[4]) > 0 then 1 else 0)
   + (if min(x[5], x[3]) > 1 then 1 else 0)
   + (if max(x[2], x[3]) == 1 then 1 else 0)
   >= 3;
