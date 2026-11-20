var x{1..4} integer >= 0 <= 2;

minimize o:
   0*x[1] - 4*x[2] + x[3] + 6*x[4] + 3*(if max(x[2], x[2]) > 2 then 1 else 0);

s.t. c1: (max(x[2], x[2]) > 2) ==> ((max(x[2], x[2]) > 2) or (x[2] >= x[3]));
s.t. c2: (max(x[2], x[2]) > 2) ==> (x[2] >= x[3]);
s.t. c3:
   (if max(x[2], x[2]) > 2 then 1 else 0)
   + (if max(x[2], x[2]) > 2 then 1 else 0) + (if x[2] >= x[3] then 1 else 0)
   >= 3;
