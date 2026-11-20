var x{1..6} integer >= 0 <= 1;
maximize o: (-4)*x[1] + (-3)*x[2] + (4)*x[3] + (3)*x[4] + (5)*x[5] + (2)*x[6];
s.t. c1: ((if x[4] != 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[4], x[3], x[6]);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
