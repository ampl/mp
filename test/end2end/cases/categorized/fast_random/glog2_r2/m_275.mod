var x{1..5} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (4)*x[2] + (4)*x[3] + (3)*x[4] + (-3)*x[5];
s.t. c1: not ((x[5] == 4) and (x[3] >= 4));
s.t. c2: ((if x[2] < 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((if x[3] >= 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c4: (x[2] >= 1) ==> (x[2] < 1) else (x[2] >= 4);
s.t. c5: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) > 0;
s.t. c6: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
