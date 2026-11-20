var x{1..6} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (3)*x[2] + (5)*x[3] + (2)*x[4] + (6)*x[5] + (-5)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= 2);
