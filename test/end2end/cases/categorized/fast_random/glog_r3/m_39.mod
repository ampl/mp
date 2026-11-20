var x{1..6} integer >= -2 <= 0;
minimize o: (4)*x[1] + (1)*x[2] + (6)*x[3] + (-5)*x[4] + (4)*x[5] + (-1)*x[6];
s.t. c1: (x[5] != -1) and (x[4] != -1) ==> (x[4] >= -1);
s.t. c2: ((if x[1] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
