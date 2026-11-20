var x{1..4} integer >= 0 <= 3;
maximize o: (2)*x[1] + (1)*x[2] + (4)*x[3] + (-1)*x[4] + (-3)*(if x[4] <= 3 then 1 else 0);
s.t. c1: (x[4] <= 3) and ((x[2] < x[2]) or not (abs(x[3] - x[1]) >= 3));
s.t. c2: (if abs(x[3] - x[1]) >= 3 then 1 else 0) + (if abs(x[3] - x[1]) >= 3 then 1 else 0) + (if abs(x[3] - x[1]) >= 3 then 1 else 0) >= 1;
