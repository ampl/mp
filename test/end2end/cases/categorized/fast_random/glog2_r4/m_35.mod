var x{1..6} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (2)*x[2] + (-1)*x[3] + (-3)*x[4] + (1)*x[5] + (-4)*x[6];
s.t. c1: (x[6] < 0) ==> (x[4] <= 1) else (x[2] <= 1);
s.t. c2: (x[3] < 2) ==> (x[5] < 2) else (x[3] == 0);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 0;
s.t. c4: alldiff(x[2], x[3], x[4], x[6]);
s.t. c5: (x[5] != 2) and (x[1] <= 2) ==> (x[4] >= 2);
s.t. c6: exists {i in 1..6} x[i] == 2;
