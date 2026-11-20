var x{1..4} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (1)*x[2] + (-5)*x[3] + (5)*x[4];
s.t. c1: (x[3] != 0) ==> (x[3] != 1);
s.t. c2: (x[4] == 0) ==> (x[4] > 0);
s.t. c3: ((if x[2] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
