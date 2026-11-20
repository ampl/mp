var x{1..6} integer >= -2 <= 0;
minimize o: (-3)*x[1] + (1)*x[2] + (-4)*x[3] + (-4)*x[4] + (-5)*x[5] + (3)*x[6];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 0;
s.t. c2: (x[6] <= 0) ==> (x[4] > -2) else (x[1] >= 0);
s.t. c3: (x[3] > -1) or (x[4] != -2) or (x[5] != -2);
s.t. c4: (x[6] >= 0) or (x[4] == -1) or (x[4] == 0);
s.t. c5: ((if x[3] == -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
