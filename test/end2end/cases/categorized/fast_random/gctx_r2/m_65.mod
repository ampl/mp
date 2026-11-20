var x{1..4} integer >= -1 <= 1;
maximize o: (-1)*x[1] + (5)*x[2] + (3)*x[3] + (0)*x[4] + (-4)*(if abs(x[4] - x[2]) == 2 then 1 else 0);
s.t. c1: (if min(x[4], x[1]) < -1 then x[3] else x[4]) < 1;
s.t. c2: (abs(x[4] - x[2]) == 2) and ((abs(x[4] - x[2]) == 2) or not (min(x[4], x[1]) < -1));
