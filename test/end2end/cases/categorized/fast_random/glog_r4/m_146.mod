var x{1..3} integer >= 0 <= 2;
minimize o: (5)*x[1] + (5)*x[2] + (-1)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: (x[3] < 1) <==> (x[1] + x[1] <= 1);
