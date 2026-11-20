var x{1..4} integer >= -1 <= 1;
maximize o: (4)*x[1] + (-3)*x[2] + (4)*x[3] + (0)*x[4];
s.t. c1: (x[1] != x[3]) and ((x[1] + x[2] == 1) or not (x[4] > -1));
s.t. c2: (x[1] + x[2] == 1) or (x[1] != x[3]);
