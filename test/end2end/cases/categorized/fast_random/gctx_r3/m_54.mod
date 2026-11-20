var x{1..3} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-4)*x[2] + (0)*x[3];
s.t. c1: not (max(x[2], x[1]) < 3);
s.t. c2: not (x[2] + x[3] > 1);
s.t. c3: not ((abs(x[2] - x[3]) <= 3) and (x[2] == 0));
s.t. c4: (abs(x[2] - x[3]) <= 3) or (max(x[2], x[1]) < 3);
