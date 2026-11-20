var x{1..5} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (-1)*x[2] + (-4)*x[3] + (1)*x[4] + (1)*x[5];
s.t. c1: max(x[5], x[5]) - min(x[2], x[1]) < 1;
s.t. c2: alldiff(x[3], x[1]);
