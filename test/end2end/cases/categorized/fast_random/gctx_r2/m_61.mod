var x{1..4} integer >= -1 <= 1;

minimize o:
   -x[1] + 0*x[2] + 6*x[3] + 6*x[4] + 4*(if min(x[4], x[2]) == 1 then 1 else 0);

s.t. c1:
   (if max(x[1], x[2]) != -1 then 1 else 0)
   + (if max(x[1], x[2]) != -1 then 1 else 0)
   + (if min(x[4], x[2]) == 1 then 1 else 0)
   >= 3;
s.t. c2:
   (min(x[4], x[2]) == 1)
   ==> ((min(x[4], x[2]) == 1) or (min(x[4], x[2]) == 1));
s.t. c3:
   (min(x[4], x[2]) == 1)
   ==> ((min(x[4], x[2]) == 1) or (max(x[1], x[2]) != -1));
