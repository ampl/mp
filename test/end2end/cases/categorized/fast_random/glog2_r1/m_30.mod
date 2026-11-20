var x{1..4} integer >= 0 <= 4;
maximize o: (4)*x[1] + (-3)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c2: max(x[4], x[4]) - min(x[1], x[3]) <= 4;
s.t. c3: alldiff(x[1], x[2], x[3]);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
