var x{1..3} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-1)*x[2] + (-1)*x[3];
s.t. c1: (x[2] > 2) or (x[3] != 2) or (x[1] != 1);
