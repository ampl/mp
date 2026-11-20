var x{1..3} integer >= -1 <= 2;
maximize o: (5)*x[1] + (-1)*x[2] + (-5)*x[3] + (1)*(if x[1] <= 0 then 1 else 0);
s.t. c1: (x[1] <= 0) <==> (abs(x[1] - x[1]) == 3);
s.t. c2: (if x[2] >= x[3] then 1 else 0) + (if x[2] >= x[3] then 1 else 0) + (if x[2] >= x[3] then 1 else 0) >= 1;
