var x{1..3} integer >= -1 <= 1;
minimize o: (-4)*x[1] + (2)*x[2] + (1)*x[3] + (-6)*(if abs(x[2] - x[2]) == 0 then 1 else 0);
s.t. c1: (if abs(x[2] - x[2]) == 0 then x[1] else 0) + (if abs(x[2] - x[2]) == 0 then x[1] else 0) >= -1;
s.t. c2: (abs(x[2] - x[2]) == 0) or (abs(x[2] - x[2]) == 0);
s.t. c3: (if abs(x[2] - x[1]) <= 2 then x[2] else x[2]) > 0;
s.t. c4: (abs(x[2] - x[2]) == 0) <==> (abs(x[2] - x[1]) <= 2);
