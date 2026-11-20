var x{1..3} integer >= 0 <= 3;

minimize o: x[1] - x[2] - x[3];

s.t. c1: not ((x[1] == 3) and (x[1] == 3));
s.t. c2: (x[1] == 3) ==> (x[1] == 3);
s.t. c3:
   (min(x[3], x[1]) != 1)
   and ((abs(x[1] - x[2]) != 1) or not (abs(x[3] - x[2]) >= 2));
s.t. c4:
   (if min(x[3], x[1]) != 1 then x[2] else 0)
   + (if abs(x[1] - x[2]) != 1 then x[3] else 0)
   > 0;
