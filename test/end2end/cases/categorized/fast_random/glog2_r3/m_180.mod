var x{1..3} integer >= 0 <= 4;
maximize o: (2)*x[1] + (-1)*x[2] + (-4)*x[3];
s.t. c1: ((if x[2] < 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c2: (x[3] >= 1) ==> (x[3] == 0) else (x[2] <= 3);
