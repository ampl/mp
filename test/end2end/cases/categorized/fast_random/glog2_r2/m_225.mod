var x{1..4} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (0)*x[2] + (2)*x[3] + (-4)*x[4];
s.t. c1: (x[4] != 0) ==> (x[4] < 3);
s.t. c2: ((if x[2] != 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c3: (x[4] < 3) ==> (x[3] < 1) else (x[1] >= 2);
s.t. c4: (x[2] != 3) or (x[1] <= 1) or (x[3] != 0);
s.t. c5: (if x[4] > 0 then x[4] else x[3]) <= 2;
