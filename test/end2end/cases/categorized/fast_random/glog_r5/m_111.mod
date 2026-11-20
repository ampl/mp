var x{1..6} integer >= 0 <= 4;
maximize o: (2)*x[1] + (4)*x[2] + (5)*x[3] + (4)*x[4] + (2)*x[5] + (1)*x[6];
s.t. c1: alldiff(x[1], x[6], x[2]);
