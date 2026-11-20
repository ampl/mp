var x{1..3} integer >= -2 <= 1;
maximize o: (-666)*x[1] + (2000)*x[2] + (-999)*x[3];
s.t. c1: max(x[3], x[1]) - min(x[3], x[3]) > -2;
s.t. c2: abs(x[1] - x[1]) == 1;
