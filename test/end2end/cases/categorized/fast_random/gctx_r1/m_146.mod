var x{1..5} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (3)*x[2] + (5)*x[3] + (-5)*x[4] + (1)*x[5];
s.t. c1: not ((x[2] + x[1] >= 2) and (x[3] <= 1));
s.t. c2: (x[3] <= 1) ==> (x[3] <= 1);
s.t. c3: (if x[2] + x[1] >= 2 then 1 else 0) + (if x[3] <= 1 then 1 else 0) + (if x[3] <= 1 then 1 else 0) >= 3;
