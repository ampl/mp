var x{1..4} integer >= -2 <= 1;
maximize o: (-1)*x[1] + (5)*x[2] + (-4)*x[3] + (-4)*x[4];
s.t. c1: (if x[2] < 0 then x[4] else x[3]) < 0;
