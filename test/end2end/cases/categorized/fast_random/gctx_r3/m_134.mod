var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] - 3*x[2] + x[3];

s.t. c1: (x[3] > 0) and ((x[3] > 0) or not (x[3] > 0));
s.t. c2:
   (if min(x[3], x[3]) == 2 then x[3] else 0)
   + (if min(x[3], x[3]) == 2 then x[2] else 0)
   >= 2;
s.t. c3:
   (min(x[3], x[3]) == 2)
   ==> ((min(x[3], x[3]) == 2) or (min(x[3], x[3]) == 2));
s.t. c4: (min(x[3], x[3]) == 2) ==> ((x[3] > 0) or (x[3] > 0));
