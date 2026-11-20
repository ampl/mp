var x{1..5} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (-4)*x[2] + (-1)*x[3] + (-5)*x[4] + (4)*x[5];
s.t. c1: not ((x[2] <= 4) and (x[1] != 1));
