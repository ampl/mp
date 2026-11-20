var x{1..5} integer >= -1 <= 2;
minimize o: (-3)*x[1] + (-5)*x[2] + (2)*x[3] + (1)*x[4] + (2)*x[5];
s.t. c1: (if x[4] + x[1] > 0 then 1 else 0) + (if abs(x[5] - x[5]) >= 3 then 1 else 0) + (if abs(x[5] - x[5]) >= 3 then 1 else 0) >= 3;
s.t. c2: not (x[4] + x[1] > 0);
