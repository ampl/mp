var x{1..5} integer >= -3 <= 1;
maximize o: (-4)*x[1] + (-4)*x[2] + (2)*x[3] + (2)*x[4] + (-3)*x[5];
s.t. c1: abs(x[3] - x[2]) < 2;
