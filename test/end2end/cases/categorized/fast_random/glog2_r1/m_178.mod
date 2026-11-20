var x{1..7} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (6)*x[2] + (-2)*x[3] + (3)*x[4] + (2)*x[5] + (0)*x[6] + (3)*x[7];
s.t. c1: (if x[4] < 0 then x[7] else x[5]) < 0;
s.t. c2: abs(x[7] - x[5]) <= 1;
s.t. c3: (if x[3] <= 0 then x[7] else x[4]) >= 1;
