var x{1..4} integer >= 0 <= 2;
maximize o: (4)*x[1] + (5)*x[2] + (-1)*x[3] + (-3)*x[4] + (0)*(if x[4] > 2 then 1 else 0);
s.t. c1: not (abs(x[4] - x[4]) < 0);
s.t. c2: (x[4] > 2) ==> (x[4] > 2);
s.t. c3: not ((x[1] < 1) and (x[1] < 1));
s.t. c4: (x[4] > 2) <==> (x[1] < 1);
