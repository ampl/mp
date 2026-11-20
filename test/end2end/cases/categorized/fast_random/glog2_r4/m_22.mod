var x{1..4} integer >= 0 <= 2;
minimize o: (3)*x[1] + (5)*x[2] + (-2)*x[3] + (-2)*x[4];
s.t. c1: (x[1] <= 2) ==> (x[4] < 2) else (x[4] == 1);
s.t. c2: (if x[1] < 1 then x[1] else x[2]) > 2;
