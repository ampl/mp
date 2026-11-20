var x{1..3} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (6)*x[2] + (-3)*x[3];
s.t. c1: alldiff(x[3], x[1]);
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
