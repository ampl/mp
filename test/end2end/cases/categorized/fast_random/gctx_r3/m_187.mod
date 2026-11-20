var x{1..5} integer >= 0 <= 3;
minimize o: (5)*x[1] + (0)*x[2] + (-3)*x[3] + (-2)*x[4] + (-2)*x[5];
s.t. c1: (if abs(x[4] - x[4]) <= 2 then x[5] else 0) + (if abs(x[4] - x[4]) <= 2 then x[4] else 0) <= 3;
s.t. c2: not ((abs(x[4] - x[4]) <= 2) and (min(x[4], x[4]) >= 0));
s.t. c3: (min(x[4], x[4]) >= 0) ==> ((min(x[4], x[4]) >= 0) or (min(x[4], x[4]) >= 0));
