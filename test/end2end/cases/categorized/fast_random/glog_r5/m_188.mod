var x{1..3} integer >= -2 <= 2;
maximize o: (-3)*x[1] + (-3)*x[2] + (-4)*x[3];
s.t. c1: (x[1] == -2) and (x[2] >= 1) ==> (x[2] >= -2);
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: (x[2] >= -1) ==> (x[3] != -1) else (x[2] >= 0);
