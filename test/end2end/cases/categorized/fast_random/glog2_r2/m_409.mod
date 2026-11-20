var x{1..4} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-1)*x[2] + (-4)*x[3] + (5)*x[4];
s.t. c1: abs(x[4] - x[1]) >= 1;
s.t. c2: exists {i in 1..4} x[i] == 2;
s.t. c3: abs(x[4] - x[1]) >= 1;
s.t. c4: exists {i in 1..4} x[i] == 2;
