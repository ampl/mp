var x{1..5} integer >= 0 <= 3;
minimize o: (6)*x[1] + (6)*x[2] + (4)*x[3] + (-3)*x[4] + (-4)*x[5];
s.t. c1: abs(x[5] - x[3]) < 0;
s.t. c2: (if x[2] < 1 then x[2] else x[2]) > 1;
