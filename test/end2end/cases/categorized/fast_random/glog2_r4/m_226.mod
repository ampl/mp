var x{1..3} integer >= 0 <= 3;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: (x[1] > 0) <==> (x[1] + x[3] > 3);
