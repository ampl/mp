var x{1..5} integer >= -2 <= 3;
minimize o: (-5)*x[1] + (-2)*x[2] + (-2)*x[3] + (6)*x[4] + (-2)*x[5];
s.t. c1: max(x[4], x[1]) - min(x[5], x[4]) <= -1;
