var x{1..4} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (4)*x[2] + (-3)*x[3] + (1)*x[4];
s.t. c1: (x[1] == -1) ==> (x[4] > -1);
s.t. c2: exists {i in 1..4} x[i] == -2;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: exists {i in 1..4} x[i] == -2;
