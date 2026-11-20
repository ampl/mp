var x{1..4} integer >= 0 <= 1;
minimize o: (-4)*x[1] + (-4)*x[2] + (2)*x[3] + (0)*x[4];
s.t. c1: (if x[3] < 0 then x[1] else x[2]) > 1;
s.t. c2: not ((x[2] < 1) and (x[3] <= 0));
