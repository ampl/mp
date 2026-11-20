var x{1..3} integer >= -2 <= 1;
maximize o: (3)*x[1] + (0)*x[2] + (-5)*x[3];
s.t. c1: (x[1] != -2) and (x[1] > -1) ==> (x[1] >= -2);
s.t. c2: ((if x[2] == -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3]) == 3;
s.t. c4: abs(x[3] - x[1]) == 3;
