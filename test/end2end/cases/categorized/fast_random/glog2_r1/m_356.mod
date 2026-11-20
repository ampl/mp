var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c2: (x[4] != 1) ==> (x[1] != 2);
s.t. c3: (x[3] <= 0) <==> (x[1] + x[1] == 2);
s.t. c4: ((if x[3] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
