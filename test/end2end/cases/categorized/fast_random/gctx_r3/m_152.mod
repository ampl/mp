var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-2)*x[2] + (-3)*x[3] + (3)*(if x[2] < 1 then 1 else 0);
s.t. c1: (abs(x[2] - x[3]) < 0) ==> (min(x[1], x[2]) == 2) else (abs(x[2] - x[3]) < 0);
s.t. c2: (if x[2] < 1 then 1 else 0) + (if x[2] < 1 then 1 else 0) + (if x[2] < 1 then 1 else 0) >= 1;
s.t. c3: (x[2] + x[1] != 0) and ((min(x[1], x[2]) == 2) or not (x[2] < 1));
s.t. c4: (if x[2] < 1 then x[3] else 0) + (if x[2] < 1 then x[3] else 0) >= 0;
