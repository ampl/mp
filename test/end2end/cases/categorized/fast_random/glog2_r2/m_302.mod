var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-5)*x[2] + (-4)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[3], x[3]) == 0;
