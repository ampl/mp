var x{1..5} integer >= -3 <= 1;
maximize o: (4)*x[1] + (-1)*x[2] + (-2)*x[3] + (-4)*x[4] + (6)*x[5];
s.t. c1: alldiff(x[1], x[5]);
