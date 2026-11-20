var x{1..3} integer >= 0 <= 2;
maximize o: (3)*x[1] + (-5)*x[2] + (1)*x[3] + (1)*(if x[3] < 2 then 1 else 0);
s.t. c1: (min(x[2], x[3]) == 2) ==> (x[3] < 2) else (x[3] < 0);
s.t. c2: (x[3] < 2) <==> (x[3] < 0);
s.t. c3: (if x[3] < 0 then x[3] else x[2]) > 1;
s.t. c4: not (max(x[1], x[3]) >= 2);
