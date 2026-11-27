var x{1..4} integer >= 0 <= 2;

minimize o: 4*x[1] - x[2] + 4*x[3] - 5*x[4];

s.t. c1: (x[4] >= 2) and ((max(x[2], x[4]) == 2) or not (x[4] >= 2));
s.t. c2:
   (if x[4] + x[2] <= 1 then x[1] else 0)
   + (if x[4] + x[2] <= 1 then x[3] else 0)
   >= 1;
s.t. c3:
   (if x[4] + x[2] <= 1 then x[1] else 0)
   + (if max(x[2], x[4]) == 2 then x[2] else 0)
   <= 1;
