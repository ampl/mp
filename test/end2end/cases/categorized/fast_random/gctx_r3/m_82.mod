var x{1..4} integer >= 0 <= 3;

minimize o: -5*x[1] + x[2] - 5*x[3] - 3*x[4];

s.t. c1: (if min(x[4], x[2]) > 0 then x[4] else x[4]) < 3;
s.t. c2: not (min(x[1], x[3]) == 0);
s.t. c3: not ((min(x[4], x[2]) > 0) and (min(x[1], x[3]) == 0));
s.t. c4:
   (min(x[4], x[2]) > 0) ==> ((min(x[4], x[2]) > 0) or (min(x[4], x[2]) > 0));
