var x{1..6} integer >= 0 <= 1;
maximize o: (2)*x[1] + (-1)*x[2] + (-3)*x[3] + (3)*x[4] + (-3)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[1], x[2], x[5]);
