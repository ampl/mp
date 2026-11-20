var x{1..4} integer >= -2 <= 0;
maximize o: (-3)*x[1] + (1)*x[2] + (0)*x[3] + (2)*x[4];
s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
