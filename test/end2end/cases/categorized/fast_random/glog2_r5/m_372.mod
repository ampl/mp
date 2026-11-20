var x{1..3} integer >= -2 <= 0;
maximize o: (3)*x[1] + (4)*x[2] + (1)*x[3];
s.t. c1: alldiff(x[2], x[1]);
