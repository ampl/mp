var x{1..4} integer >= -1 <= 1;

minimize o: -2*x[1] + 2*x[2] + 3*x[3] + 5*x[4];

s.t. c1: (max(x[1], x[4]) == -1) ==> (max(x[3], x[2]) <= 1);
s.t. c2:
   (if max(x[3], x[2]) <= 1 then x[3] else 0)
   + (if max(x[3], x[2]) <= 1 then x[3] else 0)
   < 0;
