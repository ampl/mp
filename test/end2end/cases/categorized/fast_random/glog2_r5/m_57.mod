var x{1..4} integer >= 0 <= 4;
minimize o: (-333)*x[1] + (-1000)*x[2] + (2000)*x[3] + (51)*x[4];
s.t. c1: not ((x[1] >= 3) and (x[4] <= 2));
s.t. c2: ((if x[4] != 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c3: atmost 1 {i in 1..4} (x[i] >= 4);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) >= 3;
