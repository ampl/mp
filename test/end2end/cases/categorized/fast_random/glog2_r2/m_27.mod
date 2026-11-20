var x{1..3} integer >= 0 <= 2;
minimize o: (-333)*x[1] + (-1000)*x[2] + (-666)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: not ((x[1] > 1) and (x[2] <= 2));
s.t. c3: (x[1] > 1) or (x[2] >= 0) or (x[2] >= 2);
