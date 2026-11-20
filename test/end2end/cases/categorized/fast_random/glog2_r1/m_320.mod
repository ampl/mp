var x{1..5} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (-1)*x[2] + (1)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: abs(x[1] - x[2]) > 3;
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 0);
s.t. c3: (if x[1] <= 0 then x[4] else x[2]) > 0;
