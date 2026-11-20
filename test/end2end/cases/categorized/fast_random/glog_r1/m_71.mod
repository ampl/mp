var x{1..6} integer >= 0 <= 4;
minimize o: (3)*x[1] + (2)*x[2] + (-4)*x[3] + (4)*x[4] + (1)*x[5] + (1)*x[6];
s.t. c1: alldiff(x[3], x[5], x[6]);
s.t. c2: exists {i in 1..6} x[i] == 4;
s.t. c3: exists {i in 1..6} x[i] == 4;
s.t. c4: (x[3] != 0) ==> (x[4] != 0);
