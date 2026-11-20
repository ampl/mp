var x{1..5} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (-4)*x[2] + (0)*x[3] + (0)*x[4] + (3)*x[5];
s.t. c1: (max(x[3], x[5]) < 0) ==> ((max(x[3], x[5]) < 0) or (x[2] >= 1));
s.t. c2: (max(x[3], x[5]) < 0) and ((x[2] >= 1) or not (x[2] >= 1));
s.t. c3: (if x[4] <= 2 then 1 else 0) + (if x[2] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0) >= 2;
s.t. c4: (if x[2] >= 1 then x[2] else x[4]) <= 2;
