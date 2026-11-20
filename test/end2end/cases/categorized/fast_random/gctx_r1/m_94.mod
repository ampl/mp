var x{1..4} integer >= -1 <= 2;
minimize o: (1)*x[1] + (3)*x[2] + (2)*x[3] + (-3)*x[4];
s.t. c1: (min(x[4], x[1]) == 1) ==> (abs(x[2] - x[4]) < 2) else (min(x[4], x[1]) == 1);
s.t. c2: (min(x[2], x[2]) <= 1) and ((min(x[2], x[2]) <= 1) or not (min(x[2], x[2]) <= 1));
s.t. c3: (min(x[4], x[1]) == 1) or (min(x[2], x[2]) <= 1);
