var x{1..5} integer >= 0 <= 4;
minimize o: (6)*x[1] + (6)*x[2] + (0)*x[3] + (0)*x[4] + (-1)*x[5];
s.t. c1: abs(x[1] - x[4]) != 0;
s.t. c2: (if x[3] == 4 then x[3] else x[3]) != 0;
s.t. c3: not ((x[1] <= 3) and (x[4] >= 0));
s.t. c4: alldiff(x[5], x[4], x[2], x[1]);
