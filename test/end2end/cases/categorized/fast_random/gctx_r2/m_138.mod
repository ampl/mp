var x{1..5} integer >= -1 <= 1;
minimize o: (6)*x[1] + (-5)*x[2] + (-5)*x[3] + (-1)*x[4] + (6)*x[5] + (-2)*(if x[2] + x[5] < 0 then 1 else 0);
s.t. c1: (if x[2] + x[5] < 0 then x[1] else x[4]) > 0;
s.t. c2: (if x[2] + x[5] < 0 then 1 else 0) + (if x[3] + x[3] >= 0 then 1 else 0) + (if x[3] + x[3] >= 0 then 1 else 0) >= 2;
