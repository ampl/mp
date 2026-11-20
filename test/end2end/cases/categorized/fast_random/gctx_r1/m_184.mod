var x{1..3} integer >= -1 <= 2;
minimize o: (4)*x[1] + (3)*x[2] + (-5)*x[3] + (2)*(if abs(x[3] - x[3]) != 2 then 1 else 0);
s.t. c1: (if x[3] < 2 then 1 else 0) + (if abs(x[3] - x[3]) != 2 then 1 else 0) + (if abs(x[3] - x[3]) != 2 then 1 else 0) >= 2;
s.t. c2: (max(x[3], x[3]) < -1) or (max(x[3], x[3]) < -1);
s.t. c3: (max(x[3], x[3]) < -1) and ((x[3] < 2) or not (abs(x[3] - x[3]) != 2));
s.t. c4: (if x[3] < 2 then x[3] else x[2]) == 1;
