var x{1..4} integer >= -3 <= 2;
minimize o: (6)*x[1] + (-1)*x[2] + (-5)*x[3] + (3)*x[4];
s.t. c1: atleast 1 {i in 1..4} (x[i] >= 1);
s.t. c2: not ((x[2] > -1) and (x[3] == -2));
