var x{1..4} integer >= 0 <= 3;
maximize o: (2)*x[1] + (3)*x[2] + (3)*x[3] + (-3)*x[4];
s.t. c1: (x[2] <= 3) or (x[1] > 3) or (x[3] == 0);
