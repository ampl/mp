var x{1..6} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (5)*x[2] + (4)*x[3] + (6)*x[4] + (-3)*x[5] + (-5)*x[6];
s.t. c1: atleast 2 {i in 1..6} (x[i] >= 1);
s.t. c2: exists {i in 1..6} x[i] == 4;
s.t. c3: x[6] == 3;
