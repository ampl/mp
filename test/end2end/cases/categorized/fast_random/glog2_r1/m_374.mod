var x{1..5} integer >= -3 <= 1;
minimize o: (-2)*x[1] + (5)*x[2] + (2)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: ((if x[1] != -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c2: not ((x[1] <= 1) and (x[5] <= 1));
s.t. c3: (x[4] >= 0) ==> (x[2] != -1);
s.t. c4: (x[4] == -2) and (x[5] >= 1) ==> (x[3] >= -2);
