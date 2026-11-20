var x{1..4} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (-5)*x[2] + (3)*x[3] + (-2)*x[4];
s.t. c1: abs(x[3] - x[1]) != 1;
