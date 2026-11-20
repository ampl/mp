var x{1..5} integer >= 0 <= 2;
maximize o: (-5)*x[1] + (0)*x[2] + (3)*x[3] + (5)*x[4] + (-4)*x[5];
s.t. c1: abs(x[1] - x[3]) == 2;
