var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] - 2*x[2] - 5*x[3];

s.t. c1: not (abs(x[2] - x[1]) >= 1);
s.t. c2: (max(x[2], x[1]) == 3) ==> ((x[3] == 3) or (min(x[3], x[3]) == 1));
s.t. c3:
   (if abs(x[2] - x[1]) >= 1 then x[3] else 0)
   + (if min(x[3], x[3]) == 1 then x[2] else 0)
   <= 1;
