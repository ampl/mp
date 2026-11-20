var x{1..4} integer >= 0 <= 3;
minimize o: (-333)*x[1] + (17)*x[2] + (-3000)*x[3] + (-666)*x[4];
s.t. c1: (x[1] > 0) or (x[1] < 0) or (x[2] > 0);
s.t. c2: max(x[1], x[1]) - min(x[2], x[1]) == 3;
s.t. c3: not ((x[1] == 3) and (x[1] <= 1));
