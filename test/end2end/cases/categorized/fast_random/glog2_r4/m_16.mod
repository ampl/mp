var x{1..4} integer >= -2 <= 0;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: abs(x[3] - x[3]) != 2;
s.t. c2: (if x[4] != 0 then x[4] else x[1]) <= -2;
s.t. c3: not ((x[4] == -1) and (x[1] != 0));
s.t. c4: (x[3] != 0) and (x[4] >= 0) ==> (x[1] >= 0);
