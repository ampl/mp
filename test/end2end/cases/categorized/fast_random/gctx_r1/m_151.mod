var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-4)*x[2] + (2)*x[3] + (1)*x[4];
s.t. c1: (if x[1] < 2 then 1 else 0) + (if abs(x[3] - x[2]) < 0 then 1 else 0) + (if abs(x[3] - x[2]) < 0 then 1 else 0) >= 2;
s.t. c2: not ((x[1] < 2) and (abs(x[3] - x[2]) < 0));
