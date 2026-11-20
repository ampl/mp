var x{1..3} integer >= 0 <= 4;
minimize o: (5)*x[1] + (-4)*x[2] + (3)*x[3];
s.t. c1: (if x[2] == 4 then x[3] else x[3]) == 3;
