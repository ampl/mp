var x{1..6} integer >= 0 <= 2;
minimize o: (-3000)*x[1] + (-333)*x[2] + (-666)*x[3] + (-666)*x[4] + (-333)*x[5] + (3000)*x[6];
s.t. c1: (x[2] < 1) ==> (x[1] > 2);
s.t. c2: ((if x[4] != 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
s.t. c4: (x[1] != 1) or (x[5] <= 2) or (x[3] < 2);
