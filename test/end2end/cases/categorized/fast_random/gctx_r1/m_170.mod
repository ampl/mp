var x{1..3} integer >= 0 <= 2;

minimize o:
   -5*x[1] - 4*x[2] - 5*x[3] + 3*(if min(x[2], x[3]) < 1 then 1 else 0);

s.t. c1:
   (if max(x[3], x[2]) >= 0 then x[3] else 0)
   + (if abs(x[1] - x[1]) < 2 then x[3] else 0)
   > 2;
s.t. c2: not ((max(x[3], x[2]) >= 0) and (x[1] != 1));
s.t. c3: (abs(x[1] - x[1]) < 2) ==> (min(x[2], x[3]) < 1);
s.t. c4: not (max(x[3], x[2]) >= 0);
