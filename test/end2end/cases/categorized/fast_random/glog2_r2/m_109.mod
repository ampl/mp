var x{1..6} integer >= -2 <= 1;
minimize o: (4)*x[1] + (-2)*x[2] + (-5)*x[3] + (0)*x[4] + (2)*x[5] + (4)*x[6];
s.t. c1: abs(x[6] - x[2]) != 2;
s.t. c2: abs(x[2] - x[3]) != 1;
s.t. c3: abs(x[6] - x[2]) != 2;
s.t. c4: abs(x[2] - x[3]) != 1;
