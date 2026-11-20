var x{1..4} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (2)*x[2] + (3)*x[3] + (-5)*x[4];
s.t. c1: max(x[2], x[2]) - min(x[4], x[2]) >= 2;
s.t. c2: exists {i in 1..4} x[i] == 0;
