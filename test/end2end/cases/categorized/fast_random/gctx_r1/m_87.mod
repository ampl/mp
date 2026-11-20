var x{1..5} integer >= -1 <= 1;
minimize o: (-2)*x[1] + (1)*x[2] + (6)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (abs(x[5] - x[5]) >= 0) <==> (abs(x[5] - x[5]) >= 0);
s.t. c2: not ((abs(x[2] - x[3]) < 0) and (x[5] + x[3] <= 0));
s.t. c3: (if abs(x[5] - x[5]) >= 0 then 1 else 0) + (if x[5] + x[3] <= 0 then 1 else 0) + (if x[5] + x[3] <= 0 then 1 else 0) >= 3;
s.t. c4: (if x[5] + x[3] <= 0 then x[4] else x[3]) < 0;
