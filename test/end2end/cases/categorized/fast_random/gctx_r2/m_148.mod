var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (1)*x[2] + (-5)*x[3];
s.t. c1: (if abs(x[3] - x[3]) != 0 then x[1] else 0) + (if abs(x[3] - x[3]) != 0 then x[2] else 0) > 2;
s.t. c2: (abs(x[3] - x[3]) != 0) or (min(x[3], x[2]) != 0);
