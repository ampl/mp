var x{1..4} integer >= -1 <= 1;
minimize o: (-2)*x[1] + (-5)*x[2] + (4)*x[3] + (5)*x[4] + (4)*(if x[3] > x[2] then 1 else 0);
s.t. c1: (if x[3] > x[2] then 1 else 0) + (if max(x[4], x[2]) > 1 then 1 else 0) + (if x[3] + x[1] > 0 then 1 else 0) >= 1;
s.t. c2: (if max(x[4], x[2]) > 1 then x[4] else x[1]) == 0;
s.t. c3: (if x[3] > x[2] then 1 else 0) + (if x[3] + x[1] > 0 then 1 else 0) + (if max(x[4], x[2]) > 1 then 1 else 0) >= 2;
s.t. c4: (max(x[4], x[2]) > 1) <==> (x[3] > x[2]);
