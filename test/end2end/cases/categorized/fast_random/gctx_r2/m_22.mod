var x{1..3} integer >= -1 <= 1;

minimize o: 6*x[1] + 3*x[2] + 5*x[3];

s.t. c1: (min(x[1], x[2]) >= 1) or (x[2] < 1);
s.t. c2:
   (if abs(x[2] - x[1]) < 2 then 1 else 0) + (if x[2] + x[1] > 1 then 1 else 0)
   + (if x[2] < 1 then 1 else 0)
   >= 1;
