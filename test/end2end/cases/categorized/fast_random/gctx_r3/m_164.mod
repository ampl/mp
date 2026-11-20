var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-1)*x[2] + (-3)*x[3] + (6)*x[4];
s.t. c1: (abs(x[4] - x[1]) <= 3) or (abs(x[4] - x[1]) <= 3);
s.t. c2: not (x[1] > 3);
s.t. c3: (abs(x[4] - x[1]) <= 3) ==> ((x[1] > 3) or (x[1] > 3));
s.t. c4: (if x[4] > x[1] then x[1] else 0) + (if x[4] > x[1] then x[3] else 0) > 3;
