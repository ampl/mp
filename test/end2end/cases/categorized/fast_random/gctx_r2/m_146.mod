var x{1..3} integer >= 0 <= 3;

minimize o: x[1] + 3*x[2] + 6*x[3];

s.t. c1:
   (if x[1] + x[3] > 2 then 1 else 0) + (if max(x[3], x[2]) != 2 then 1 else 0)
   + (if abs(x[1] - x[3]) < 1 then 1 else 0)
   >= 1;
s.t. c2: not (x[1] + x[3] > 2);
s.t. c3: not (max(x[3], x[2]) != 2);
