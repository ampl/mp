var x{1..5} integer >= 0 <= 2;
maximize o: (3)*x[1] + (6)*x[2] + (-4)*x[3] + (-4)*x[4] + (4)*x[5] + (4)*(if abs(x[2] - x[4]) < 2 then 1 else 0);
s.t. c1: (if x[2] >= x[2] then x[5] else x[2]) > 0;
s.t. c2: not (x[2] >= x[2]);
s.t. c3: (if x[2] >= x[2] then 1 else 0) + (if abs(x[2] - x[4]) < 2 then 1 else 0) + (if x[2] >= x[2] then 1 else 0) >= 2;
s.t. c4: (x[2] >= x[2]) and ((x[5] < x[5]) or not (x[2] >= x[2]));
