var x{1..4} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (-4)*x[2] + (-5)*x[3] + (4)*x[4];
s.t. c1: (x[3] == 0) ==> (x[3] > 1) else (x[4] > 0);
s.t. c2: (x[4] > 1) ==> (x[4] <= 0);
s.t. c3: (x[3] == 0) ==> (x[4] >= 1);
s.t. c4: ((if x[2] > 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
