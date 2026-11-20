var x{1..4} integer >= 0 <= 3;
minimize o: (4)*x[1] + (-1)*x[2] + (-4)*x[3] + (1)*x[4];
s.t. c1: not ((x[1] > 0) and (x[3] < 1));
s.t. c2: (x[1] > 2) ==> (x[1] >= 1) else (x[4] <= 0);
s.t. c3: alldiff(x[3], x[4]);
s.t. c4: (x[3] < 0) or (x[4] < 0) or (x[2] == 2);
s.t. c5: ((if x[2] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c6: (x[4] == 2) ==> (x[4] == 2) else (x[4] <= 1);
