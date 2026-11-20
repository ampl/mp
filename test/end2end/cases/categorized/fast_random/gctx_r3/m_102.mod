var x{1..4} integer >= 0 <= 3;
minimize o: (2)*x[1] + (-5)*x[2] + (4)*x[3] + (-4)*x[4];
s.t. c1: (if x[1] != 3 then 1 else 0) + (if x[2] != 0 then 1 else 0) + (if abs(x[2] - x[4]) >= 3 then 1 else 0) >= 2;
s.t. c2: (abs(x[2] - x[4]) >= 3) and ((x[1] != 3) or not (abs(x[2] - x[4]) >= 3));
s.t. c3: (x[2] != 0) ==> (x[2] != 0) else (x[1] != 3);
