var x{1..6} integer >= 0 <= 3;
minimize o: (4)*x[1] + (-3)*x[2] + (4)*x[3] + (3)*x[4] + (-1)*x[5] + (-3)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= 3);
s.t. c2: abs(x[1] - x[3]) < 2;
s.t. c3: exists {i in 1..6} x[i] == 1;
