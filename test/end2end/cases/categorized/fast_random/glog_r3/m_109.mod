var x{1..5} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (4)*x[2] + (1)*x[3] + (-2)*x[4] + (-4)*x[5];
s.t. c1: max(x[2], x[4]) - min(x[4], x[1]) == 1;
