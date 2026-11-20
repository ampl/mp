var x{1..5} integer >= -2 <= 0;
minimize o: (-4)*x[1] + (-4)*x[2] + (-1)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: abs(x[1] - x[5]) <= 0;
