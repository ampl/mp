var x{1..5} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (0)*x[2] + (-5)*x[3] + (1)*x[4] + (5)*x[5];
s.t. c1: not (abs(x[2] - x[4]) == 0);
s.t. c2: (x[5] < 2) or (abs(x[2] - x[4]) == 0);
s.t. c3: (abs(x[2] - x[4]) == 0) and ((abs(x[5] - x[4]) < 3) or not (x[5] < 2));
