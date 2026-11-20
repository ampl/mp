var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] - 3*x[2] + 4*x[3] + (if x[2] >= 3 then 1 else 0);

s.t. c1: not ((x[1] > x[1]) and (x[2] >= 3));
s.t. c2:
   (if x[3] == x[2] then x[3] else 0) + (if x[3] == x[2] then x[3] else 0) >= 1;
