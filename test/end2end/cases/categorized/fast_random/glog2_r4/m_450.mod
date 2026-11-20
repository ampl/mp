var x{1..4} integer >= -2 <= 1;
minimize o: (5)*x[1] + (-3)*x[2] + (-4)*x[3] + (-2)*x[4];
s.t. c1: alldiff(x[3], x[2]);
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 0);
s.t. c4: (x[2] != 0) ==> (x[3] <= -2) else (x[3] >= 0);
