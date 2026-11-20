var x{1..7} integer >= 0 <= 1;
minimize o: (-4)*x[1] + (6)*x[2] + (-5)*x[3] + (3)*x[4] + (5)*x[5] + (-4)*x[6] + (-5)*x[7];
s.t. c1: (if x[1] <= 1 then x[1] else x[5]) < 1;
s.t. c2: exists {i in 1..7} x[i] == 0;
s.t. c3: (if x[5] != 1 then x[1] else x[1]) < 0;
