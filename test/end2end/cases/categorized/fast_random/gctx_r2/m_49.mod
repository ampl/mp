var x{1..4} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (0)*x[2] + (2)*x[3] + (3)*x[4] + (-3)*(if x[1] + x[2] >= 0 then 1 else 0);
s.t. c1: not ((x[1] + x[2] >= 0) and (x[1] + x[2] >= 0));
s.t. c2: (abs(x[4] - x[2]) > 1) and ((abs(x[4] - x[2]) > 1) or not (abs(x[4] - x[2]) > 1));
s.t. c3: (if abs(x[4] - x[2]) > 1 then x[3] else 0) + (if abs(x[4] - x[2]) > 1 then x[3] else 0) > 1;
