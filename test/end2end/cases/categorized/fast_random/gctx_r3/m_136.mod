var x{1..5} integer >= -1 <= 2;

maximize o:
   2*x[1] + 5*x[2] + x[3] + 0*x[4] + 2*x[5]
   + 3*(if abs(x[5] - x[4]) > 0 then 1 else 0);

s.t. c1: (abs(x[5] - x[4]) > 0) ==> ((x[2] >= x[4]) or (min(x[3], x[5]) > 2));
s.t. c2:
   (if min(x[3], x[5]) > 2 then 1 else 0)
   + (if abs(x[5] - x[4]) > 0 then 1 else 0) + (if x[2] >= x[4] then 1 else 0)
   >= 1;
s.t. c3: not ((abs(x[1] - x[1]) <= 3) and (min(x[3], x[5]) > 2));
s.t. c4: not ((abs(x[5] - x[4]) > 0) and (abs(x[1] - x[1]) <= 3));
