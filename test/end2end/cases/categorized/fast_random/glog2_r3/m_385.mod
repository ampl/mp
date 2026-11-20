var x{1..4} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (4)*x[2] + (5)*x[3] + (6)*x[4];
s.t. c1: not ((x[2] < 0) and (x[3] <= 2));
s.t. c2: (x[1] > 0) <==> (x[3] + x[3] > 0);
s.t. c3: max(x[2], x[1]) - min(x[2], x[2]) <= 1;
s.t. c4: ((if x[4] != 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c5: atmost 3 {i in 1..4} (x[i] >= 2);
