var x{1..5} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-4)*x[2] + (-5)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: (x[3] > x[3]) ==> ((abs(x[1] - x[4]) != 0) or (abs(x[1] - x[4]) != 0));
s.t. c2: (abs(x[1] - x[4]) != 0) or (abs(x[1] - x[4]) != 0);
s.t. c3: (if max(x[2], x[3]) >= 1 then x[5] else 0) + (if x[4] + x[2] > 0 then x[4] else 0) < 2;
s.t. c4: not ((abs(x[1] - x[4]) != 0) and (x[4] + x[2] > 0));
