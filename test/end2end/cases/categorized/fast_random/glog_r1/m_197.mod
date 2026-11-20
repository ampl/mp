var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (3)*x[2] + (-3)*x[3] + (4)*x[4] + (1)*x[5];
s.t. c1: not ((x[1] >= 2) and (x[1] <= 0));
