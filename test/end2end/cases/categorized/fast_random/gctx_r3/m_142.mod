var x{1..5} integer >= 0 <= 3;

minimize o:
   4*x[1] - 4*x[2] + 6*x[3] + 6*x[4] - 4*x[5]
   - 5*(if x[2] + x[1] != 1 then 1 else 0);

s.t. c1:
   (if x[2] + x[1] != 1 then x[1] else 0)
   + (if max(x[3], x[4]) < 3 then x[3] else 0)
   > 0;
s.t. c2: (x[2] >= x[2]) ==> (max(x[3], x[4]) < 3);
