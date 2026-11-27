var x{1..3} integer >= -1 <= 1;

minimize o: x[1] - x[2] - 4*x[3];

s.t. c1: (if max(x[1], x[2]) > 1 then x[3] else x[1]) > 1;
s.t. c2:
   (max(x[1], x[1]) == 0)
   and ((max(x[1], x[2]) > 1) or not (max(x[1], x[2]) != 0));
s.t. c3: (max(x[1], x[1]) == 0) ==> (max(x[1], x[1]) == 0);
s.t. c4: (max(x[1], x[2]) > 1) <==> (max(x[1], x[2]) > 1);
