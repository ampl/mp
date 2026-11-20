var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (5)*x[2] + (0)*x[3];
s.t. c1: (if max(x[1], x[2]) >= 1 then x[2] else 0) + (if max(x[3], x[1]) != 1 then x[3] else 0) <= 1;
s.t. c2: (max(x[3], x[1]) != 1) <==> (max(x[3], x[1]) != 1);
s.t. c3: (max(x[1], x[2]) >= 1) <==> (max(x[1], x[2]) >= 1);
