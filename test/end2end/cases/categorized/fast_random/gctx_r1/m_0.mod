var x{1..4} integer >= -1 <= 1;

minimize o: -x[1] - 5*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[4] > 1) ==> ((x[4] > 1) or (abs(x[4] - x[3]) <= 2));
s.t. c2:
   (if abs(x[4] - x[3]) <= 2 then x[2] else 0)
   + (if abs(x[4] - x[3]) <= 2 then x[2] else 0)
   <= -1;
