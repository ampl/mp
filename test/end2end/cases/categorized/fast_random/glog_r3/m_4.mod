var x{1..4} integer >= 0 <= 4;
maximize o: (2)*x[1] + (-5)*x[2] + (-1)*x[3] + (2)*x[4];
s.t. c1: not ((x[2] >= 0) and (x[4] == 3));
s.t. c2: not ((x[4] >= 1) and (x[4] != 0));
