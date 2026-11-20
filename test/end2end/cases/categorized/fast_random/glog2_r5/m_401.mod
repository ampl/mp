var x{1..5} integer >= 0 <= 3;
maximize o: (3)*x[1] + (-4)*x[2] + (0)*x[3] + (-3)*x[4] + (1)*x[5];
s.t. c1: alldiff(x[1], x[3], x[5], x[2]);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (if x[5] < 3 then x[3] else x[1]) <= -1;
