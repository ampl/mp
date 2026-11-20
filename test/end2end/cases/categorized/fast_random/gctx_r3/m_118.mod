var x{1..4} integer >= 0 <= 2;

minimize o:
   -3*x[1] - 3*x[2] - 5*x[3] + 4*x[4]
   - 6*(if abs(x[1] - x[4]) >= 1 then 1 else 0);

s.t. c1: (min(x[4], x[2]) == 1) <==> (x[2] == 2);
s.t. c2:
   (if min(x[4], x[2]) == 1 then x[4] else 0)
   + (if abs(x[1] - x[4]) >= 1 then x[1] else 0)
   < 0;
s.t. c3: not ((x[2] == 2) and (x[3] + x[3] != 2));
