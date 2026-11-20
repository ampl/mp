var x{1..4} integer >= 0 <= 2;
minimize o: (2)*x[1] + (1)*x[2] + (5)*x[3] + (-1)*x[4];
s.t. c1: (x[1] <= 0) ==> (x[3] == 1) else (x[3] < 1);
