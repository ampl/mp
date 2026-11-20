var x{1..5} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-4)*x[2] + (1)*x[3] + (1)*x[4] + (2)*x[5];
s.t. c1: (abs(x[4] - x[4]) <= 1) or (abs(x[4] - x[4]) <= 1);
s.t. c2: not ((abs(x[4] - x[4]) <= 1) and (x[3] == x[3]));
