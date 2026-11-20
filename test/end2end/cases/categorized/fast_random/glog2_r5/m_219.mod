var x{1..5} integer >= 0 <= 2;
minimize o: (0)*x[1] + (2)*x[2] + (-3)*x[3] + (-4)*x[4] + (-5)*x[5];
s.t. c1: ((if x[2] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[4] > 2) or (x[2] == 0) or (x[1] >= 1);
s.t. c3: ((if x[2] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c4: (x[4] > 2) or (x[2] == 0) or (x[1] >= 1);
