var x{1..3} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (5)*x[2] + (-1)*x[3];
s.t. c1: (x[1] <= 3) or (x[1] <= 3);
s.t. c2: (x[1] + x[1] == 0) <==> (x[1] <= 3);
s.t. c3: (if x[1] <= 3 then 1 else 0) + (if x[1] + x[1] == 0 then 1 else 0) + (if x[1] + x[1] == 0 then 1 else 0) >= 1;
