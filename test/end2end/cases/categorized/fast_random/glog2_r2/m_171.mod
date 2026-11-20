var x{1..5} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-3)*x[2] + (-4)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: abs(x[4] - x[1]) <= 0;
