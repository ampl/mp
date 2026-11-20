var x{1..3} integer >= 0 <= 3;
maximize o: (1)*x[1] + (-5)*x[2] + (-5)*x[3];
s.t. c1: (x[2] >= 2) <==> (x[2] >= 2);
s.t. c2: (x[2] >= 2) ==> ((min(x[1], x[3]) >= 2) or (x[2] >= 2));
s.t. c3: (if x[2] >= 2 then 1 else 0) + (if x[3] == 2 then 1 else 0) + (if min(x[1], x[3]) >= 2 then 1 else 0) >= 3;
