var x{1..5} integer >= 0 <= 4;
maximize o: (5)*x[1] + (1)*x[2] + (1)*x[3] + (-1)*x[4] + (-4)*x[5];
s.t. c1: ((if x[3] == -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c3: ((if x[1] != 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
s.t. c4: (x[5] > -1) ==> (x[5] > -1) else (x[5] != 4);
