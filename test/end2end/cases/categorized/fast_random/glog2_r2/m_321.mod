var x{1..5} integer >= 0 <= 2;
minimize o: (4)*x[1] + (3)*x[2] + (0)*x[3] + (3)*x[4] + (1)*x[5];
s.t. c1: not ((x[2] <= 0) and (x[4] <= 3));
s.t. c2: not ((x[5] != 2) and (x[1] > 2));
s.t. c3: (if x[4] <= 3 then x[2] else x[3]) == 0;
s.t. c4: not ((x[1] <= 2) and (x[1] <= 0));
