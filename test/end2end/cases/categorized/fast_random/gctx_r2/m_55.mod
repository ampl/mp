var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (4)*x[2] + (-1)*x[3] + (-3)*(if abs(x[3] - x[2]) <= 0 then 1 else 0);
s.t. c1: (abs(x[3] - x[2]) <= 0) <==> (abs(x[3] - x[2]) <= 0);
s.t. c2: (max(x[1], x[3]) < 2) and ((max(x[1], x[3]) < 2) or not (max(x[1], x[3]) < 2));
s.t. c3: (if abs(x[3] - x[2]) <= 0 then x[2] else 0) + (if abs(x[3] - x[2]) <= 0 then x[1] else 0) < 0;
s.t. c4: (min(x[2], x[2]) < 2) ==> ((min(x[2], x[2]) < 2) or (min(x[2], x[2]) < 2));
