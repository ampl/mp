var x{1..4} integer >= -3 <= 1;
maximize o: (-5)*x[1] + (5)*x[2] + (-5)*x[3] + (2)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -3;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 0);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: alldiff(x[1], x[2], x[3], x[4]);
