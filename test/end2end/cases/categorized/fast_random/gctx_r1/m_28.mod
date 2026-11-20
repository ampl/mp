var x{1..4} integer >= 0 <= 3;

minimize o:
   x[1] - 3*x[2] + 4*x[3] - 4*x[4] + 2*(if max(x[3], x[2]) < 3 then 1 else 0);

s.t. c1:
   (if max(x[3], x[2]) < 3 then 1 else 0)
   + (if max(x[3], x[2]) < 3 then 1 else 0)
   + (if abs(x[3] - x[2]) >= 3 then 1 else 0)
   >= 2;
s.t. c2: (if abs(x[3] - x[2]) >= 3 then x[1] else x[1]) >= 1;
s.t. c3:
   (x[4] <= x[3]) ==> (abs(x[3] - x[2]) >= 3) else (abs(x[3] - x[2]) >= 3);
s.t. c4: (if max(x[3], x[2]) < 3 then x[1] else x[3]) <= 0;
