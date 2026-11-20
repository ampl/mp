var x{1..4} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (-4)*x[2] + (-5)*x[3] + (6)*x[4];
s.t. c1: (x[1] < 1) ==> (x[1] != 1) else (x[2] != -2);
s.t. c2: (x[2] < -2) ==> (x[4] == -2);
s.t. c3: ((if x[3] > 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c4: max(x[2], x[3]) - min(x[3], x[4]) == -1;
s.t. c5: (x[1] < 1) ==> (x[1] != 1) else (x[2] != -2);
s.t. c6: max(x[2], x[3]) - min(x[3], x[4]) == -1;
