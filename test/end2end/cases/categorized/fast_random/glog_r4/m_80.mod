var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (-5)*x[2] + (0)*x[3];
s.t. c1: ((x[3] >= 0) ==> (x[2] != 2)) and ((x[2] != 2) ==> (x[3] >= 0));
s.t. c2: ((if x[1] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: ((x[1] < 1) ==> (x[2] != 2)) and ((x[2] != 2) ==> (x[2] >= 1));
s.t. c4: numberof 2 in (x[1], x[2], x[3]) != 1;
