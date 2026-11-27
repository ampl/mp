var x{1..3} integer >= -1 <= 2;

minimize o: -x[1] + 0*x[2] + x[3] - (if max(x[1], x[3]) > 2 then 1 else 0);

s.t. c1: (max(x[1], x[3]) > 2) and ((x[2] <= 0) or not (max(x[1], x[3]) > 2));
s.t. c2: not (x[2] <= 0);
s.t. c3: (x[2] <= 0) ==> ((max(x[1], x[3]) > 2) or (max(x[1], x[3]) > 2));
