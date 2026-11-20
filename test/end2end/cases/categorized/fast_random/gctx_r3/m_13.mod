var x{1..3} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (1)*x[2] + (0)*x[3] + (3)*(if abs(x[1] - x[3]) == 0 then 1 else 0);
s.t. c1: (if abs(x[3] - x[1]) > 2 then x[3] else x[1]) <= 2;
s.t. c2: not ((abs(x[1] - x[3]) == 0) and (abs(x[1] - x[3]) == 0));
s.t. c3: (min(x[2], x[2]) > 1) ==> ((abs(x[1] - x[3]) == 0) or (abs(x[1] - x[3]) == 0));
s.t. c4: not ((min(x[2], x[2]) > 1) and (abs(x[3] - x[1]) > 2));
