var x{1..6} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: not ((x[1] < 3) and (x[3] > 4));
s.t. c2: (if x[6] >= 2 then x[3] else x[1]) > 0;
s.t. c3: abs(x[2] - x[2]) >= 1;
