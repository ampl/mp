var x{1..3} integer >= -2 <= 0;
minimize o: (6)*x[1] + (4)*x[2] + (-5)*x[3];
s.t. c1: ((if x[1] < 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] <= 0) and (x[2] > -2) ==> (x[2] >= 0);
