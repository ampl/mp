var x{1..4} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (5)*x[2] + (3)*x[3] + (-2)*x[4];
s.t. c1: abs(x[3] - x[3]) >= 1;
s.t. c2: (if x[3] == 2 then x[3] else x[1]) >= 1;
