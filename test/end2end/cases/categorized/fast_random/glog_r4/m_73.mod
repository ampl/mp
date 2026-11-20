var x{1..4} integer >= 0 <= 2;
maximize o: (2)*x[1] + (0)*x[2] + (0)*x[3] + (2)*x[4];
s.t. c1: (x[2] >= 1) or (x[2] != 0) or (x[3] == 2);
