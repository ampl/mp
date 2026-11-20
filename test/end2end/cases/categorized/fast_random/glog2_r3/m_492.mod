var x{1..5} integer >= 0 <= 4;
maximize o: (-2)*x[1] + (-3)*x[2] + (-1)*x[3] + (5)*x[4] + (4)*x[5];
s.t. c1: not ((x[5] <= 0) and (x[4] < 1));
