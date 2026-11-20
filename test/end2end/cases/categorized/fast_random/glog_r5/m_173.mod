var x{1..4} integer >= -2 <= 2;
maximize o: (5)*x[1] + (-1)*x[2] + (2)*x[3] + (3)*x[4];
s.t. c1: (if x[1] < 2 then x[4] else x[4]) != 2;
