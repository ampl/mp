var x{1..3} integer >= 0 <= 3;
minimize o: (2)*x[1] + (6)*x[2] + (-2)*x[3] + (-1)*(if abs(x[2] - x[1]) < 2 then 1 else 0);
s.t. c1: (abs(x[2] - x[1]) < 2) or (abs(x[2] - x[1]) <= 1);
s.t. c2: not (abs(x[2] - x[1]) <= 1);
s.t. c3: (abs(x[2] - x[2]) != 3) ==> (abs(x[2] - x[1]) < 2) else (abs(x[2] - x[1]) < 2);
