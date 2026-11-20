var x{1..4} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c2: (x[1] < -1) <==> (x[2] + x[2] < 0);
