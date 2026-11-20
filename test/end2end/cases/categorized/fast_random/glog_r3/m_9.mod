var x{1..6} integer >= -2 <= 0;
minimize o: (1)*x[1] + (6)*x[2] + (-3)*x[3] + (4)*x[4] + (0)*x[5] + (2)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= -1);
s.t. c2: not ((x[5] == -1) and (x[3] == 0));
s.t. c3: abs(x[4] - x[5]) != 1;
s.t. c4: not ((x[4] < 0) and (x[6] == 0));
