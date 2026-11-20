var x{1..5} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (4)*x[2] + (5)*x[3] + (-2)*x[4] + (3)*x[5] + (4)*(if abs(x[2] - x[5]) > 2 then 1 else 0);
s.t. c1: (x[1] != 3) ==> (min(x[2], x[4]) < 2) else (x[1] != 3);
s.t. c2: (abs(x[2] - x[5]) > 2) or (x[2] + x[3] != 2);
