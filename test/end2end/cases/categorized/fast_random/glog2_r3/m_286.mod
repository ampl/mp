var x{1..4} integer >= -2 <= 1;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[1], x[2], x[3], x[4]);
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 0);
s.t. c3: (x[1] == -1) ==> (x[2] == 1) else (x[2] == 1);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) <= 1;
