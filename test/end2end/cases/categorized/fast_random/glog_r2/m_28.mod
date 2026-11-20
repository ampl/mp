var x{1..4} integer >= 0 <= 4;
maximize o: (5)*x[1] + (4)*x[2] + (1)*x[3] + (-2)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[2] <= 1) ==> (x[3] < 0);
s.t. c3: (if x[3] > 4 then x[1] else x[2]) > 4;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
