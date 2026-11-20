var x{1..4} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] - 5*(if x[3] == 2 then 1 else 0);

s.t. c1: (x[2] > x[4]) ==> (x[2] > x[4]);
s.t. c2: (if x[3] == 2 then x[4] else x[1]) >= 1;
s.t. c3:
   (if abs(x[4] - x[3]) > 2 then 1 else 0)
   + (if abs(x[4] - x[3]) > 2 then 1 else 0)
   + (if abs(x[4] - x[3]) > 2 then 1 else 0)
   >= 1;
s.t. c4:
   (if x[2] > x[4] then 1 else 0) + (if x[2] > x[4] then 1 else 0)
   + (if x[2] > x[4] then 1 else 0)
   >= 2;
