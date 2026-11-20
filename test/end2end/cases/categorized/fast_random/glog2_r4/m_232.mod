var x{1..5} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (-3)*x[2] + (-4)*x[3] + (-3)*x[4] + (-2)*x[5];
s.t. c1: abs(x[5] - x[3]) < 4;
