var x{1..6} integer >= -2 <= 1;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: (x[5] > 1) and (x[3] != 1) ==> (x[3] >= 1);
s.t. c2: ((if x[4] == -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c3: ((if x[3] <= -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c4: (x[2] <= 1) ==> (x[4] == -1) else (x[5] > -1);
