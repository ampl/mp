var x{1..5} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-2)*x[2] + (-2)*x[3] + (-3)*x[4] + (-3)*x[5];
s.t. c1: (if x[1] <= 2 then x[2] else x[1]) < 2;
s.t. c2: not ((min(x[5], x[5]) < 1) and (abs(x[5] - x[4]) != 1));
s.t. c3: (if min(x[5], x[5]) < 1 then x[2] else x[3]) < 0;
