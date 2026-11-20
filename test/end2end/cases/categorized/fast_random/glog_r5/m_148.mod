var x{1..6} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (6)*x[2] + (-5)*x[3] + (5)*x[4] + (0)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: atleast 3 {i in 1..6} (x[i] >= 1);
s.t. c3: abs(x[5] - x[6]) <= 2;
