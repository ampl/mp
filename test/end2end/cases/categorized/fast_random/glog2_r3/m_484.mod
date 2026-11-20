var x{1..5} integer >= -3 <= 1;
minimize o: (6)*x[1] + (3)*x[2] + (4)*x[3] + (0)*x[4] + (-1)*x[5];
s.t. c1: atleast 2 {i in 1..5} (x[i] >= -1);
