var x{1..5} integer >= -3 <= 2;
maximize o: (4)*x[1] + (4)*x[2] + (1)*x[3] + (-4)*x[4] + (4)*x[5];
s.t. c1: abs(x[5] - x[2]) < 0;
