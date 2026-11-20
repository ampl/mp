var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-1)*x[2] + (4)*x[3] + (5)*x[4];
s.t. c1: (x[4] == 1) ==> (x[4] == 2) else (x[3] < 3);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: abs(x[3] - x[4]) >= 0;
