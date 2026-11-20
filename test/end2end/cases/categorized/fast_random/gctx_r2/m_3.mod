var x{1..5} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (-4)*x[2] + (3)*x[3] + (-4)*x[4] + (1)*x[5] + (-5)*(if x[3] + x[2] < 1 then 1 else 0);
s.t. c1: (x[5] <= 1) ==> ((x[3] + x[2] < 1) or (x[5] <= 1));
s.t. c2: (if x[5] <= 1 then 1 else 0) + (if x[3] + x[2] < 1 then 1 else 0) + (if x[5] <= 1 then 1 else 0) >= 1;
s.t. c3: not (x[5] <= 1);
