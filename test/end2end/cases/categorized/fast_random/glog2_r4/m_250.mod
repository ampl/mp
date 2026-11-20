var x{1..4} integer >= 0 <= 3;
maximize o: (3)*x[1] + (-3)*x[2] + (-4)*x[3] + (2)*x[4];
s.t. c1: (if x[3] == 2 then x[1] else x[4]) == 1;
