var x{1..4} integer >= -2 <= 1;
maximize o: (-5)*x[1] + (-2)*x[2] + (4)*x[3] + (-3)*x[4];
s.t. c1: (x[2] <= -2) or (x[2] <= 0) or (x[2] <= 1);
s.t. c2: ((if x[2] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[3] == 0) ==> (x[3] < -1);
