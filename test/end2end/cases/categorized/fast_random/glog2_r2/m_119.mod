var x{1..5} integer >= -2 <= 1;
maximize o: (4)*x[1] + (5)*x[2] + (-4)*x[3] + (5)*x[4] + (-1)*x[5];
s.t. c1: abs(x[3] - x[2]) != 1;
s.t. c2: abs(x[3] - x[2]) != 1;
