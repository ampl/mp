var x{1..6} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[2], x[3], x[4]);
