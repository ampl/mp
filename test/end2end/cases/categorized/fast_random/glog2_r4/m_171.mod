var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (2)*x[2] + (0)*x[3] + (-3)*x[4];
s.t. c1: (if x[1] > 4 then x[4] else x[1]) != 0;
