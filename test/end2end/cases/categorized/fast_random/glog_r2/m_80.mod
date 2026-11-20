var x{1..3} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (-4)*x[2] + (2)*x[3];
s.t. c1: (if x[2] <= 0 then x[1] else x[1]) > 1;
