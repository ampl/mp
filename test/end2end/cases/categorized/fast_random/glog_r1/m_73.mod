var x{1..4} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (-2)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: alldiff(x[4], x[3], x[2]);
s.t. c2: (x[2] > 2) or (x[1] > 0) or (x[3] < 1);
s.t. c3: exists {i in 1..4} x[i] == 3;
