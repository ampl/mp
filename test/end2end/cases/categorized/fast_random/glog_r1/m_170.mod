var x{1..3} integer >= -2 <= 0;
maximize o: (5)*x[1] + (3)*x[2] + (-4)*x[3];
s.t. c1: abs(x[3] - x[1]) <= 1;
