var x{1..4} integer >= -2 <= 2;
minimize o: (2)*x[1] + (-4)*x[2] + (-3)*x[3] + (0)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c2: (if x[3] != -2 then x[1] else x[2]) != -1;
s.t. c3: (x[2] == 2) ==> (x[4] != 2) else (x[1] > -1);
s.t. c4: ((if x[4] <= -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c5: atmost 3 {i in 1..4} (x[i] >= 0);
s.t. c6: (x[4] != 2) or (x[3] != 2) or (x[1] > 0);
s.t. c7: numberof 0 in (x[1], x[2], x[3], x[4]) > 1;
