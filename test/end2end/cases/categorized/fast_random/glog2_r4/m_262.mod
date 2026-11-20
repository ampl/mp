var x{1..4} integer >= 0 <= 4;
minimize o: (3)*x[1] + (2)*x[2] + (5)*x[3] + (6)*x[4];
s.t. c1: ((if x[1] > 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: not ((x[4] != 0) and (x[4] == 1));
s.t. c3: alldiff(x[2], x[3], x[1], x[4]);
s.t. c4: (x[2] <= 3) and (x[4] <= 4) ==> (x[4] >= 3);
