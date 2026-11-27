var x{1..5} integer >= 0 <= 3;

minimize o: x[1] + x[2] + 0*x[3] + 3*x[4] + 5*x[5];

s.t. c1: (if x[4] != x[3] then x[5] else x[1]) <= 2;
s.t. c2:
   (if min(x[4], x[2]) == 1 then 1 else 0) + (if x[4] != x[3] then 1 else 0)
   + (if min(x[4], x[5]) != 1 then 1 else 0)
   >= 2;
