var x{1..3} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (1)*x[2] + (-2)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: alldiff(x[3], x[1], x[2]);
