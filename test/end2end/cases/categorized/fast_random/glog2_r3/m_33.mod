var x{1..3} integer >= 0 <= 3;
maximize o: (3)*x[1] + (-5)*x[2] + (-1)*x[3];
s.t. c1: ((if x[2] <= 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c2: numberof 3 in (x[1], x[2], x[3]) < 0;
s.t. c3: max(x[3], x[3]) - min(x[1], x[1]) > 0;
s.t. c4: (x[1] == 2) and (x[1] < 1) ==> (x[1] >= 2);
s.t. c5: x[2] == 0;
