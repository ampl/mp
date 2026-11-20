var x{1..4} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (x[4] > 2) ==> (x[4] <= 0) else (x[1] <= 0);
s.t. c2: max(x[1], x[3]) - min(x[2], x[2]) > 1;
