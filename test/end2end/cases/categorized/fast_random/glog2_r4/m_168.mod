var x{1..4} integer >= 0 <= 1;
maximize o: (-4)*x[1] + (-1)*x[2] + (-3)*x[3] + (5)*x[4];
s.t. c1: (x[1] != 0) or (x[2] == 1) or (x[2] < 0);
s.t. c2: (if x[1] <= 1 then x[4] else x[2]) < 0;
