var x{1..4} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (if x[4] > -1 then x[1] else x[1]) != 2;
s.t. c2: (x[1] < 1) or (x[4] >= 1) or (x[2] != 0);
