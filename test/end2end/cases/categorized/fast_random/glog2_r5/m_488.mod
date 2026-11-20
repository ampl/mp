var x{1..4} integer >= 0 <= 1;
maximize o: (-2)*x[1] + (-3)*x[2] + (0)*x[3] + (4)*x[4];
s.t. c1: (x[4] == 0) and (x[1] > 0) ==> (x[1] >= 0);
s.t. c2: abs(x[1] - x[3]) > 1;
s.t. c3: exists {i in 1..4} x[i] == 0;
