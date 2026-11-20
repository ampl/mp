var x{1..5} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (-4)*x[2] + (2)*x[3] + (4)*x[4] + (2)*x[5];
s.t. c1: (if x[5] != 2 then x[3] else x[3]) == 2;
