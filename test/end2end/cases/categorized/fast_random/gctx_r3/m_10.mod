var x{1..4} integer >= -1 <= 1;

minimize o: 6*x[1] + 6*x[2] + 4*x[3] - 4*x[4];

s.t. c1:
   (if x[2] != x[1] then x[2] else 0)
   + (if min(x[1], x[1]) != -1 then x[4] else 0)
   >= 0;
s.t. c2:
   (if x[3] != 0 then x[3] else 0) + (if min(x[1], x[1]) != -1 then x[3] else 0)
   < -1;
s.t. c3: (x[2] != x[1]) and ((min(x[1], x[1]) != -1) or not (x[3] != 0));
