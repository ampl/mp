var x{1..5} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (0)*x[2] + (2)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: not (max(x[4], x[2]) >= 2);
s.t. c2: (if max(x[4], x[2]) != 1 then x[2] else 0) + (if abs(x[3] - x[1]) >= 2 then x[1] else 0) <= 1;
s.t. c3: (max(x[4], x[2]) >= 2) <==> (abs(x[3] - x[1]) >= 2);
s.t. c4: (max(x[4], x[2]) >= 2) <==> (max(x[4], x[2]) != 1);
