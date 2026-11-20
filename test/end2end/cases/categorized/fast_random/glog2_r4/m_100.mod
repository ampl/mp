var x{1..3} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (-1)*x[2] + (2)*x[3];
s.t. c1: max(x[1], x[3]) - min(x[1], x[2]) != 2;
s.t. c2: not ((x[1] >= 1) and (x[3] != 1));
s.t. c3: alldiff(x[2], x[1], x[3]);
