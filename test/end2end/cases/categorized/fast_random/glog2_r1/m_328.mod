var x{1..7} integer >= 0 <= 1;
maximize o: (-3)*x[1] + (-1)*x[2] + (-2)*x[3] + (-3)*x[4] + (4)*x[5] + (0)*x[6] + (5)*x[7];
s.t. c1: alldiff(x[6], x[5]);
