var x{1..4} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: atleast 3 {i in 1..4} (x[i] >= 3);
s.t. c2: abs(x[1] - x[1]) < 1;
