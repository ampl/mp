var x{1..3} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-5)*x[2] + (2)*x[3] + (2)*(if x[3] <= 0 then 1 else 0);
s.t. c1: (x[3] <= 0) and ((x[3] <= 0) or not (x[3] <= 0));
s.t. c2: (x[3] <= 0) and ((max(x[2], x[3]) > 1) or not (x[3] <= 0));
s.t. c3: (x[3] <= 0) or (x[3] <= 0);
s.t. c4: (max(x[2], x[3]) > 1) ==> (x[3] <= 0) else (x[3] <= 0);
