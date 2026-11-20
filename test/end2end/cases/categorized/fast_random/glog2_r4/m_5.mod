var x{1..6} integer >= 0 <= 2;
minimize o: (6)*x[1] + (2)*x[2] + (-2)*x[3] + (3)*x[4] + (6)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[3], x[6], x[5]);
s.t. c2: abs(x[4] - x[5]) != 1;
s.t. c3: (x[5] > 2) or (x[6] > 1) or (x[4] >= 0);
s.t. c4: (x[3] == 0) or (x[4] != 2) or (x[4] >= 2);
s.t. c5: ((if x[2] != 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
