var x{1..5} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (-5)*x[2] + (-2)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: (min(x[3], x[3]) < 2) ==> (min(x[3], x[3]) < 2) else (min(x[3], x[3]) < 2);
s.t. c2: (if x[2] > x[4] then x[5] else x[4]) < 2;
s.t. c3: (min(x[3], x[3]) < 2) and ((min(x[3], x[3]) < 2) or not (x[2] > x[4]));
