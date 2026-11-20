var x{1..6} integer >= 0 <= 2;
maximize o: (4)*x[1] + (0)*x[2] + (-5)*x[3] + (6)*x[4] + (-5)*x[5] + (1)*x[6];
s.t. c1: not ((x[6] < 2) and (x[4] < 0));
