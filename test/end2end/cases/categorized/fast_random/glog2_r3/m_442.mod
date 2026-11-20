var x{1..4} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (5)*x[2] + (1)*x[3] + (4)*x[4];
s.t. c1: (x[3] > 1) or (x[3] <= 1) or (x[2] == 2);
