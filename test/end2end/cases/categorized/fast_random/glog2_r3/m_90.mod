var x{1..4} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (6)*x[2] + (1)*x[3] + (-3)*x[4];
s.t. c1: not ((x[2] < 4) and (x[3] <= 4));
