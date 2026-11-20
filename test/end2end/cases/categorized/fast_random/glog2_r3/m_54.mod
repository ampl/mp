var x{1..3} integer >= -3 <= 2;
minimize o: (2)*x[1] + (1)*x[2] + (-3)*x[3];
s.t. c1: (x[2] == 2) <==> (x[1] + x[3] > 0);
s.t. c2: numberof -2 in (x[1], x[2], x[3]) <= 3;
s.t. c3: (x[3] <= -1) <==> (x[1] + x[1] <= 0);
s.t. c4: exists {i in 1..3} x[i] == 1;
