var x{1..4} integer >= 0 <= 3;
minimize o: (5)*x[1] + (-1)*x[2] + (2)*x[3] + (-3)*x[4];
s.t. c1: not (max(x[4], x[2]) < 0);
s.t. c2: (x[2] != x[1]) ==> (max(x[4], x[2]) < 0) else (max(x[4], x[2]) < 0);
s.t. c3: (if max(x[4], x[2]) < 0 then x[3] else 0) + (if x[2] != x[1] then x[1] else 0) >= 1;
s.t. c4: (x[2] != x[1]) ==> ((max(x[4], x[2]) < 0) or (max(x[4], x[2]) < 0));
