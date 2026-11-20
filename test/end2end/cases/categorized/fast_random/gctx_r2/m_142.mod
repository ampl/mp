var x{1..3} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (4)*x[2] + (2)*x[3];
s.t. c1: (x[1] <= 0) and ((x[1] <= 0) or not (x[2] + x[1] > 0));
s.t. c2: (if x[1] <= 0 then x[3] else 0) + (if x[1] <= 0 then x[1] else 0) >= 0;
