var x{1..5} integer >= 0 <= 1;
minimize o: (1)*x[1] + (1)*x[2] + (1)*x[3] + (1)*x[4] + (6)*x[5];
s.t. c1: abs(x[4] - x[1]) < 1;
