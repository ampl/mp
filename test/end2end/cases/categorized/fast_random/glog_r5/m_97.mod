var x{1..4} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (-2)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: (x[4] > 2) ==> (x[4] < 1) else (x[1] <= 1);
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: (x[4] <= 1) or (x[4] <= 0) or (x[4] != 3);
