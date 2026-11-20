var x{1..5} integer >= 0 <= 4;
minimize o: (3)*x[1] + (6)*x[2] + (-2)*x[3] + (3)*x[4] + (1)*x[5];
s.t. c1: (if x[5] == 3 then x[5] else x[4]) < 1;
