var x{1..4} integer >= 0 <= 4;
minimize o: (-1000)*x[1] + (-666)*x[2] + (17)*x[3] + (-2000)*x[4];
s.t. c1: alldiff(x[2], x[4], x[3]);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c3: (x[1] >= 3) or (x[4] < 3) or (x[1] >= 2);
