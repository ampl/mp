var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-1)*x[2] + (4)*x[3] + (-1)*x[4] + (0)*x[5];
s.t. c1: (min(x[1], x[5]) < 2) and ((min(x[1], x[4]) == 1) or not (min(x[1], x[5]) < 2));
s.t. c2: (min(x[1], x[5]) < 2) and ((min(x[1], x[4]) == 1) or not (x[1] == x[4]));
