var x{1..4} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (0)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: not ((x[2] != x[2]) and (x[3] > 1));
s.t. c2: (x[3] > 1) or (x[2] != x[2]);
