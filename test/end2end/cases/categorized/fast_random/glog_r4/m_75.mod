var x{1..4} integer >= 0 <= 4;
maximize o: (2)*x[1] + (5)*x[2] + (5)*x[3] + (-4)*x[4];
s.t. c1: (x[4] == 2) ==> (x[1] < 1) else (x[4] > 3);
s.t. c2: abs(x[4] - x[2]) < 3;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
