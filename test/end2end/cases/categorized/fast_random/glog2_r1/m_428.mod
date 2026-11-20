var x{1..5} integer >= 0 <= 1;
minimize o: (1)*x[1] + (3)*x[2] + (-4)*x[3] + (6)*x[4] + (-1)*x[5];
s.t. c1: (x[5] < 1) or (x[1] < 1) or (x[3] > 1);
s.t. c2: not ((x[5] > 1) and (x[2] <= 1));
s.t. c3: not ((x[3] != 0) and (x[4] == 0));
