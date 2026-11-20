var x{1..4} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (-1)*x[2] + (-4)*x[3] + (-4)*x[4];
s.t. c1: (x[2] < 1) and (x[4] >= 2) ==> (x[1] >= 1);
s.t. c2: ((if x[4] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[2] == 1) ==> (x[4] != 1) else (x[3] < 2);
