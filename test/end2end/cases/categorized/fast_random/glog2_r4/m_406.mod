var x{1..6} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: ((if x[5] <= 2 then 1 else 0) + (if x[6] >= 2 then 1 else 0)) >= 1;
s.t. c2: not ((x[4] != 3) and (x[2] == 2));
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 3;
s.t. c4: (x[2] < 0) ==> (x[3] <= 3) else (x[1] == 0);
