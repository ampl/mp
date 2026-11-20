var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (6)*x[2] + (4)*x[3] + (6)*x[4] + (-2)*x[5];
s.t. c1: not ((max(x[5], x[1]) < 1) and (x[4] >= x[3]));
s.t. c2: (x[4] >= x[3]) or (max(x[5], x[1]) < 1);
