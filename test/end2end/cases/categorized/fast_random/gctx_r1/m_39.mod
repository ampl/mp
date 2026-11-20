var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (0)*x[2] + (-5)*x[3] + (-3)*x[4];
s.t. c1: not (abs(x[3] - x[3]) <= 2);
s.t. c2: (if abs(x[3] - x[3]) <= 2 then x[2] else x[2]) < 2;
