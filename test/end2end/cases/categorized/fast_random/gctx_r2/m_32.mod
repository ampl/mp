var x{1..3} integer >= 0 <= 2;
minimize o: (5)*x[1] + (0)*x[2] + (-1)*x[3];
s.t. c1: (min(x[1], x[2]) > 1) or (min(x[1], x[2]) > 1);
s.t. c2: (if min(x[1], x[2]) > 1 then 1 else 0) + (if min(x[1], x[2]) > 1 then 1 else 0) + (if abs(x[3] - x[1]) == 1 then 1 else 0) >= 2;
s.t. c3: not (abs(x[3] - x[1]) == 1);
