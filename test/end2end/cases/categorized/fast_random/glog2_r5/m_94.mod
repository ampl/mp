var x{1..5} integer >= -3 <= 1;
maximize o: (5)*x[1] + (1)*x[2] + (-3)*x[3] + (-5)*x[4] + (-1)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[5], x[2]) >= -2;
s.t. c2: not ((x[1] > 1) and (x[2] == 0));
s.t. c3: abs(x[1] - x[2]) <= 3;
