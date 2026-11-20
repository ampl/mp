var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (1)*x[2] + (4)*x[3] + (-5)*x[4] + (-3)*x[5];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c2: ((x[1] >= 0) ==> (x[2] != 2)) and ((x[2] != 2) ==> (x[3] >= 0));
s.t. c3: ((if x[4] <= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: atleast 3 {i in 1..5} (x[i] >= 3);
