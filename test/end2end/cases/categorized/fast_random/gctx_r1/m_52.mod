var x{1..3} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (-3)*x[2] + (-2)*x[3] + (6)*(if x[1] < x[2] then 1 else 0);
s.t. c1: (x[1] < x[2]) <==> (x[1] < 1);
s.t. c2: (x[1] < 1) or (x[3] >= 0);
s.t. c3: (x[1] < x[2]) or (x[1] < x[2]);
s.t. c4: (if x[1] < x[2] then 1 else 0) + (if x[3] >= 0 then 1 else 0) + (if x[1] < x[2] then 1 else 0) >= 1;
