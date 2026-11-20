var x{1..6} integer >= -2 <= 0;
minimize o: (-1)*x[1] + (-5)*x[2] + (0)*x[3] + (6)*x[4] + (1)*x[5] + (-1)*x[6];
s.t. c1: (x[4] > -2) ==> (x[4] == -2);
s.t. c2: alldiff(x[4], x[6], x[2], x[5]);
s.t. c3: atmost 3 {i in 1..6} (x[i] >= -1);
s.t. c4: (x[4] > -2) ==> (x[4] == -2);
s.t. c5: atmost 3 {i in 1..6} (x[i] >= -1);
