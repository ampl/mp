var x{1..5} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (-1)*x[2] + (3)*x[3] + (-2)*x[4] + (5)*x[5];
s.t. c1: (if x[5] >= 0 then x[3] else x[4]) < 1;
s.t. c2: abs(x[3] - x[3]) > 2;
s.t. c3: abs(x[5] - x[4]) >= 0;
s.t. c4: not ((x[2] > 1) and (x[1] > 2));
