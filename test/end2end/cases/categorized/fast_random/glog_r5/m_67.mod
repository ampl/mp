var x{1..6} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (4)*x[2] + (-2)*x[3] + (2)*x[4] + (-2)*x[5] + (-5)*x[6];
s.t. c1: not ((x[2] == 2) and (x[4] > -2));
s.t. c2: (x[5] <= 0) and (x[5] < -1) ==> (x[5] >= 0);
s.t. c3: alldiff(x[4], x[5], x[1]);
s.t. c4: exists {i in 1..6} x[i] == 0;
