var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-4)*x[2] + (4)*x[3] + (-5)*(if x[3] <= 2 then 1 else 0);
s.t. c1: (x[1] >= 0) ==> ((x[1] >= 0) or (x[3] < 0));
s.t. c2: (if x[3] > 2 then 1 else 0) + (if x[3] < 0 then 1 else 0) + (if x[3] <= 2 then 1 else 0) >= 3;
s.t. c3: not ((x[3] < 0) and (x[3] < 0));
s.t. c4: (x[3] <= 2) ==> (x[3] > 2) else (x[3] > 2);
