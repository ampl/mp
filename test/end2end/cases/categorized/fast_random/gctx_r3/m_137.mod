var x{1..5} integer >= 0 <= 2;
maximize o: (3)*x[1] + (2)*x[2] + (-3)*x[3] + (0)*x[4] + (3)*x[5];
s.t. c1: not ((abs(x[3] - x[2]) < 0) and (x[3] >= 1));
s.t. c2: not (x[1] >= 1);
s.t. c3: (if x[3] >= 1 then x[1] else x[3]) >= 2;
s.t. c4: (x[3] >= 1) <==> (x[3] >= 1);
