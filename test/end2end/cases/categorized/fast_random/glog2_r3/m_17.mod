var x{1..4} integer >= 0 <= 4;
minimize o: (3000)*x[1] + (-999)*x[2] + (-666)*x[3] + (17)*x[4];
s.t. c1: max(x[4], x[1]) - min(x[3], x[3]) != 2;
