var x{1..6} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (2)*x[2] + (1)*x[3] + (5)*x[4] + (5)*x[5] + (5)*x[6];
s.t. c1: alldiff(x[3], x[6], x[5], x[2]);
