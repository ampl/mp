var x{1..4} integer >= -3 <= 2;
maximize o: (1)*x[1] + (1)*x[2] + (3)*x[3] + (4)*x[4];
s.t. c1: not ((x[4] > 0) and (x[4] == -3));
