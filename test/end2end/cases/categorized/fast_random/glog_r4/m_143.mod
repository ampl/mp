var x{1..4} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (6)*x[2] + (-3)*x[3] + (-3)*x[4];
s.t. c1: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c2: (x[4] > 1) ==> (x[4] < 2) else (x[4] != 0);
s.t. c3: max(x[3], x[4]) - min(x[1], x[4]) == 3;
