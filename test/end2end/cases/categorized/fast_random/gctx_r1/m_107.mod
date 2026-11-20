var x{1..4} integer >= 0 <= 2;
minimize o: (1)*x[1] + (3)*x[2] + (0)*x[3] + (-3)*x[4];
s.t. c1: (abs(x[2] - x[1]) < 0) and ((abs(x[4] - x[4]) == 0) or not (abs(x[4] - x[4]) == 0));
s.t. c2: (abs(x[2] - x[1]) < 0) and ((abs(x[4] - x[4]) == 0) or not (abs(x[2] - x[1]) < 0));
s.t. c3: (abs(x[2] - x[1]) < 0) and ((abs(x[2] - x[1]) < 0) or not (abs(x[2] - x[1]) < 0));
