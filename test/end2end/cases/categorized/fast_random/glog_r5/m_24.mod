var x{1..4} integer >= 0 <= 3;
maximize o: (4)*x[1] + (-4)*x[2] + (0)*x[3] + (4)*x[4];
s.t. c1: (if x[2] != 2 then x[2] else x[4]) == 3;
