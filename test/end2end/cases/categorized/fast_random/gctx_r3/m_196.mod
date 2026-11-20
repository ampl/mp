var x{1..5} integer >= -1 <= 2;
minimize o: (1)*x[1] + (-4)*x[2] + (0)*x[3] + (6)*x[4] + (-3)*x[5];
s.t. c1: not ((abs(x[5] - x[1]) < 1) and (x[5] > 0));
s.t. c2: (if abs(x[5] - x[1]) < 1 then x[4] else 0) + (if max(x[5], x[5]) <= 1 then x[5] else 0) <= 0;
