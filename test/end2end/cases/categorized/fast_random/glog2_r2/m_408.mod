var x{1..4} integer >= 0 <= 1;
maximize o: (-5)*x[1] + (-4)*x[2] + (1)*x[3] + (-3)*x[4];
s.t. c1: abs(x[2] - x[1]) <= 1;
