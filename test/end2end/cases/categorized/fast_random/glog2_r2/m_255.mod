var x{1..5} integer >= -2 <= 2;
minimize o: (1)*x[1] + (6)*x[2] + (-2)*x[3] + (4)*x[4] + (0)*x[5];
s.t. c1: (if x[4] >= -2 then x[1] else x[1]) <= -2;
s.t. c2: (x[3] != 0) ==> (x[5] <= 2);
s.t. c3: not ((x[3] != 1) and (x[5] < 2));
s.t. c4: alldiff(x[1], x[2], x[4], x[3]);
s.t. c5: not ((x[4] == -2) and (x[4] >= -1));
s.t. c6: (x[5] >= 0) ==> (x[1] <= -1) else (x[3] < 0);
s.t. c7: (x[2] >= -2) or (x[5] <= 0) or (x[2] < 1);
