var x{1..3} integer >= 0 <= 2;
maximize o: (-5)*x[1] + (5)*x[2] + (-2)*x[3];
s.t. c1: (if x[3] <= 0 then x[3] else 0) + (if min(x[3], x[3]) >= 2 then x[2] else 0) <= 1;
s.t. c2: (if x[3] <= 0 then x[2] else 0) + (if abs(x[1] - x[1]) < 0 then x[2] else 0) < 1;
s.t. c3: (if max(x[2], x[1]) != 0 then x[2] else x[3]) >= 2;
s.t. c4: (abs(x[1] - x[1]) < 0) and ((min(x[3], x[3]) >= 2) or not (abs(x[1] - x[1]) < 0));
