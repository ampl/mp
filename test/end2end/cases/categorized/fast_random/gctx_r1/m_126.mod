var x{1..4} integer >= 0 <= 3;

minimize o:
   -2*x[1] - 2*x[2] + x[3] + x[4] - 4*(if abs(x[2] - x[3]) >= 1 then 1 else 0);

s.t. c1:
   (if x[3] != 1 then x[3] else 0) + (if abs(x[3] - x[3]) == 1 then x[1] else 0)
   > 3;
s.t. c2: not ((abs(x[2] - x[3]) >= 1) and (x[4] <= 2));
s.t. c3:
   (if abs(x[3] - x[3]) == 1 then x[4] else 0)
   + (if abs(x[3] - x[3]) == 1 then x[2] else 0)
   <= 3;
s.t. c4:
   (if abs(x[2] - x[3]) >= 1 then 1 else 0) + (if x[4] <= 2 then 1 else 0)
   + (if abs(x[3] - x[3]) == 1 then 1 else 0)
   >= 2;
