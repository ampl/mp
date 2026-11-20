var x{1..6} integer >= -2 <= 0;
maximize o: (1)*x[1] + (2)*x[2] + (-3)*x[3] + (1)*x[4] + (2)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[3], x[5], x[4], x[1]);
