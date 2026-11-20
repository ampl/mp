var x{1..4} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (-5)*x[2] + (5)*x[3] + (6)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: max(x[2], x[2]) - min(x[2], x[2]) == -2;
s.t. c3: (if x[2] <= 2 then x[2] else x[2]) == -1;
