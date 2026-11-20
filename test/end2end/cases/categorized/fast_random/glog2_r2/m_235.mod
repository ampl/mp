var x{1..6} integer >= -2 <= 0;
minimize o: (-1)*x[1] + (3)*x[2] + (4)*x[3] + (6)*x[4] + (5)*x[5] + (3)*x[6];
s.t. c1: ((if x[4] < -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c2: (x[6] == 0) ==> (x[3] > -1) else (x[5] < -2);
s.t. c3: not ((x[3] < 0) and (x[3] == 0));
s.t. c4: (x[1] == -2) ==> (x[4] != -1);
s.t. c5: atmost 1 {i in 1..6} (x[i] >= -2);
s.t. c6: alldiff(x[4], x[3]);
