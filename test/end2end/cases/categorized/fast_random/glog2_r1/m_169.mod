var x{1..4} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (-1)*x[2] + (-5)*x[3] + (-4)*x[4];
s.t. c1: ((x[2] > 0) ==> (x[4] > 1)) and ((x[4] > 1) ==> (x[4] >= 0));
s.t. c2: ((if x[4] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((x[2] > 0) ==> (x[4] > 1)) and ((x[4] > 1) ==> (x[4] >= 0));
s.t. c4: ((if x[4] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
