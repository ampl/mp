var x{1..4} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (4)*x[2] + (-4)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[3], x[2], x[4]);
s.t. c2: not ((x[4] >= 0) and (x[3] != 2));
