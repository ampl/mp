var x{1..5} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (6)*x[2] + (-5)*x[3] + (-5)*x[4] + (6)*x[5];
s.t. c1: (x[4] > x[3]) ==> (x[3] > 2) else (x[4] > x[3]);
s.t. c2: (x[4] != x[4]) and ((x[3] > 2) or not (x[3] > 2));
s.t. c3: (if x[3] > 2 then x[5] else 0) + (if x[3] > 2 then x[1] else 0) < 2;
s.t. c4: (if x[4] != x[4] then 1 else 0) + (if x[3] > 2 then 1 else 0) + (if x[3] > 2 then 1 else 0) >= 3;
