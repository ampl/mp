var x{1..4} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (5)*x[2] + (-2)*x[3] + (-1)*x[4];
s.t. c1: alldiff(x[1], x[3]);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[1] != 2) or (x[4] != 4) or (x[3] >= 2);
