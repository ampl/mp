var x{1..4} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-2)*x[2] + (-3)*x[3] + (-4)*x[4];
s.t. c1: ((if x[4] < 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: (x[4] <= 2) and (x[4] != 3) ==> (x[2] >= 2);
s.t. c3: (x[2] == 3) ==> (x[4] >= 2) else (x[1] <= 1);
s.t. c4: x[2] == 3;
s.t. c5: x[1] == 1;
