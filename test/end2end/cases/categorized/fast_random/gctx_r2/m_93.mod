var x{1..5} integer >= 0 <= 3;
minimize o: (3)*x[1] + (-3)*x[2] + (-3)*x[3] + (-2)*x[4] + (6)*x[5];
s.t. c1: (min(x[5], x[1]) > 1) <==> (abs(x[3] - x[4]) > 1);
s.t. c2: (min(x[5], x[1]) > 1) or (abs(x[3] - x[4]) > 1);
s.t. c3: not ((max(x[2], x[1]) == 0) and (x[5] < 3));
