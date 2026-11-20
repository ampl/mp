var x{1..5} integer >= -2 <= 1;
maximize o: (0)*x[1] + (-4)*x[2] + (1)*x[3] + (-4)*x[4] + (2)*x[5];
s.t. c1: max(x[5], x[1]) - min(x[3], x[2]) != 1;
s.t. c2: alldiff(x[2], x[4], x[1], x[3]);
s.t. c3: x[5] == 0;
