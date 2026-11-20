var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (-1)*x[2] + (-1)*x[3];
s.t. c1: (if x[3] <= 2 then x[2] else 0) + (if min(x[1], x[3]) != 2 then x[2] else 0) <= 2;
s.t. c2: (x[1] <= x[3]) <==> (x[3] <= 2);
s.t. c3: not (min(x[1], x[3]) != 2);
s.t. c4: not (min(x[1], x[3]) != 2);
