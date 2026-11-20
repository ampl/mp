var x{1..6} integer >= 0 <= 3;
minimize o: (4)*x[1] + (0)*x[2] + (1)*x[3] + (1)*x[4] + (0)*x[5] + (3)*x[6];
s.t. c1: (x[1] == 1) and (x[5] <= 2) ==> (x[3] >= 1);
s.t. c2: ((if x[2] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c3: (x[5] > 3) ==> (x[4] != 2) else (x[3] > 0);
s.t. c4: not ((x[4] >= 2) and (x[1] <= 3));
