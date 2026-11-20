var x{1..4} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (0)*x[2] + (-3)*x[3] + (-3)*x[4];
s.t. c1: (x[2] <= 3) or (x[4] < 1) or (x[1] < 3);
s.t. c2: abs(x[4] - x[2]) >= 2;
s.t. c3: (x[2] <= 3) or (x[4] < 1) or (x[1] < 3);
s.t. c4: abs(x[4] - x[2]) >= 2;
