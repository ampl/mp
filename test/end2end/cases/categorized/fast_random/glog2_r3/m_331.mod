var x{1..3} integer >= 0 <= 4;
minimize o: (3)*x[1] + (1)*x[2] + (2)*x[3];
s.t. c1: (if x[1] > 4 then x[1] else x[3]) <= 0;
