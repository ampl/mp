var x{1..5} integer >= -2 <= 2;
maximize o: (-1)*x[1] + (-1)*x[2] + (-1)*x[3] + (4)*x[4] + (1)*x[5];
s.t. c1: alldiff(x[2], x[1], x[5], x[4]);
