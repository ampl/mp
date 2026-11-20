var x{1..3} integer >= 0 <= 2;
minimize o: (1)*x[1] + (3)*x[2] + (6)*x[3];
s.t. c1: not ((x[2] <= 1) and (x[2] == 1));
s.t. c2: not ((x[1] < 2) and (x[3] == 0));
