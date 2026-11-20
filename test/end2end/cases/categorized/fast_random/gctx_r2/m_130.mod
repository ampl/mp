var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (0)*x[2] + (-1)*x[3];
s.t. c1: not (min(x[3], x[3]) != 2);
s.t. c2: (if min(x[3], x[3]) != 2 then x[1] else 0) + (if x[2] >= x[1] then x[2] else 0) < 0;
