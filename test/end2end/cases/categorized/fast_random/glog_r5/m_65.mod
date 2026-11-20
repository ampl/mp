var x{1..5} integer >= 0 <= 2;
minimize o: (3)*x[1] + (-5)*x[2] + (-2)*x[3] + (6)*x[4] + (1)*x[5];
s.t. c1: (x[1] <= 2) ==> (x[2] != 0);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c3: ((if x[3] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..5} x[i] == 0;
