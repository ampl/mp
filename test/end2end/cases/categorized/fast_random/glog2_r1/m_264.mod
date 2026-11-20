var x{1..4} integer >= -2 <= 3;
maximize o: (2)*x[1] + (-1)*x[2] + (-5)*x[3] + (3)*x[4];
s.t. c1: not ((x[2] < 3) and (x[4] >= 2));
