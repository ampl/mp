var x{1..4} integer >= 0 <= 2;
minimize o: (2)*x[1] + (3)*x[2] + (0)*x[3] + (-5)*x[4];
s.t. c1: (x[3] >= 1) or (abs(x[4] - x[2]) < 0);
s.t. c2: (max(x[4], x[3]) >= 1) ==> (abs(x[4] - x[2]) < 0);
