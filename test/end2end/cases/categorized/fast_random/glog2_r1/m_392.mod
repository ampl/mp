var x{1..3} integer >= 0 <= 2;
minimize o: (5)*x[1] + (-4)*x[2] + (-2)*x[3];
s.t. c1: (x[3] >= 0) ==> (x[2] > 2) else (x[2] >= 2);
