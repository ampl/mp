var x{1..4} integer >= 0 <= 2;
maximize o: (1)*x[1] + (1)*x[2] + (-2)*x[3] + (-1)*x[4];
s.t. c1: not ((x[2] == 1) and (x[1] < 2));
