var x{1..4} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-1)*x[2] + (-2)*x[3] + (0)*x[4];
s.t. c1: (abs(x[2] - x[3]) >= 2) or (x[2] == 1);
s.t. c2: (x[1] == 0) and ((abs(x[2] - x[3]) >= 2) or not (x[2] == 1));
