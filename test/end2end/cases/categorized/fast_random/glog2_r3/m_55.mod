var x{1..4} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-1)*x[2] + (0)*x[3] + (4)*x[4];
s.t. c1: max(x[1], x[3]) - min(x[3], x[2]) < 2;
s.t. c2: (if x[3] > 1 then x[2] else x[1]) < 1;
s.t. c3: (x[3] > 0) <==> (x[2] + x[1] != 1);
s.t. c4: max(x[2], x[3]) - min(x[4], x[4]) != 2;
s.t. c5: (x[4] <= 2) ==> (x[1] < 0) else (x[1] != 1);
s.t. c6: (x[2] <= 1) and (x[4] <= 0) ==> (x[1] >= 1);
