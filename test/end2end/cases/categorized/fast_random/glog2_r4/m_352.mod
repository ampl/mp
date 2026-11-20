var x{1..4} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (1)*x[2] + (6)*x[3] + (6)*x[4];
s.t. c1: (x[3] == 3) or (x[1] == 0) or (x[4] > 0);
