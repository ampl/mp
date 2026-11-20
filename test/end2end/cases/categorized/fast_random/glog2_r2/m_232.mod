var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (5)*x[2] + (-3)*x[3] + (5)*x[4];
s.t. c1: (x[3] <= 1) ==> (x[3] >= 1) else (x[3] > 0);
s.t. c2: max(x[3], x[1]) - min(x[1], x[3]) != 0;
s.t. c3: abs(x[1] - x[2]) == 2;
