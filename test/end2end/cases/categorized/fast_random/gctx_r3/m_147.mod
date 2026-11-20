var x{1..5} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (-2)*x[2] + (-1)*x[3] + (6)*x[4] + (6)*x[5];
s.t. c1: (x[3] <= x[3]) or (x[4] + x[3] <= 2);
s.t. c2: not (x[3] <= x[3]);
s.t. c3: (if x[4] + x[3] <= 2 then x[2] else x[2]) <= 1;
