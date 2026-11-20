var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (6)*x[2] + (-3)*x[3] + (5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: abs(x[1] - x[1]) >= 1;
s.t. c4: exists {i in 1..4} x[i] == 2;
