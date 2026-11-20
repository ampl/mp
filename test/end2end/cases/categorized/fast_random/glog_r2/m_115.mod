var x{1..3} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (6)*x[2] + (-4)*x[3];
s.t. c1: max(x[1], x[3]) - min(x[1], x[2]) == 1;
s.t. c2: (x[3] != 3) or (x[3] == 1) or (x[1] >= 4);
s.t. c3: max(x[1], x[3]) - min(x[2], x[1]) > 0;
