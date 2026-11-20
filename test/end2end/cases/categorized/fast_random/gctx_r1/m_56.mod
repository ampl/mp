var x{1..5} integer >= -1 <= 2;
maximize o: (-4)*x[1] + (0)*x[2] + (-2)*x[3] + (3)*x[4] + (-3)*x[5];
s.t. c1: (if abs(x[4] - x[4]) < 3 then 1 else 0) + (if x[4] + x[1] != 1 then 1 else 0) + (if abs(x[4] - x[4]) < 3 then 1 else 0) >= 1;
s.t. c2: (x[4] + x[1] != 1) or (x[4] + x[1] != 1);
s.t. c3: (if max(x[5], x[4]) >= 0 then 1 else 0) + (if x[4] + x[1] != 1 then 1 else 0) + (if max(x[5], x[4]) >= 0 then 1 else 0) >= 1;
