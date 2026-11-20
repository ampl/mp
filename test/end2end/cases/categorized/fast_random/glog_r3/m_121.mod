var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (1)*x[2] + (-3)*x[3];
s.t. c1: (x[2] <= 0) ==> (x[2] > 0);
s.t. c2: max(x[3], x[3]) - min(x[1], x[1]) == 1;
