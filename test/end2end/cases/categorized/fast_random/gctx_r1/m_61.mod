var x{1..3} integer >= -1 <= 2;

minimize o: -3*x[1] - 3*x[2] + x[3];

s.t. c1:
   (max(x[1], x[3]) == -1) and ((x[1] == -1) or not (max(x[1], x[3]) == -1));
s.t. c2:
   (if max(x[1], x[3]) == -1 then 1 else 0)
   + (if abs(x[1] - x[3]) > 1 then 1 else 0)
   + (if max(x[1], x[3]) == -1 then 1 else 0)
   >= 1;
s.t. c3:
   (abs(x[1] - x[3]) > 1) and ((max(x[1], x[3]) == -1) or not (x[1] == -1));
s.t. c4:
   (if max(x[1], x[3]) == -1 then x[3] else 0)
   + (if x[1] == -1 then x[1] else 0)
   <= -1;
