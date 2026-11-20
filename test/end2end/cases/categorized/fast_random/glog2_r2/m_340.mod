var x{1..3} integer >= 0 <= 4;
maximize o: (4)*x[1] + (0)*x[2] + (-5)*x[3];
s.t. c1: not ((x[2] >= 1) and (x[3] < 0));
