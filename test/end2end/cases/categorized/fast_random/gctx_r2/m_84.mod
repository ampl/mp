var x{1..4} integer >= -1 <= 1;

maximize o: 4*x[1] - 4*x[2] + 0*x[3] + 3*x[4];

s.t. c1: not ((x[3] > x[4]) and (x[3] > x[4]));
s.t. c2:
   (if x[3] > x[4] then 1 else 0) + (if x[3] > x[4] then 1 else 0)
   + (if abs(x[4] - x[1]) < 1 then 1 else 0)
   >= 2;
s.t. c3: (if abs(x[4] - x[1]) < 1 then x[2] else x[3]) >= 1;
