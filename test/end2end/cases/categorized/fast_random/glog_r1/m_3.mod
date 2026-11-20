var x{1..6} integer >= 0 <= 2;
minimize o: (4)*x[1] + (4)*x[2] + (-3)*x[3] + (1)*x[4] + (1)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: abs(x[3] - x[3]) >= 1;
s.t. c3: alldiff(x[6], x[4], x[3]);
