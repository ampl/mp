var x{1..5} integer >= 0 <= 1;
maximize o: (-3)*x[1] + (1)*x[2] + (2)*x[3] + (-4)*x[4] + (6)*x[5];
s.t. c1: max(x[5], x[2]) - min(x[5], x[5]) <= 0;
