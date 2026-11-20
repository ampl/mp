var x{1..5} integer >= 0 <= 3;
minimize o: (4)*x[1] + (2)*x[2] + (1)*x[3] + (-2)*x[4] + (0)*x[5];
s.t. c1: not (max(x[4], x[4]) <= 3);
s.t. c2: (x[3] == 2) <==> (max(x[4], x[4]) <= 3);
s.t. c3: (abs(x[1] - x[5]) < 1) ==> (abs(x[1] - x[5]) < 1);
s.t. c4: (if x[3] == 2 then 1 else 0) + (if x[3] == 2 then 1 else 0) + (if abs(x[1] - x[5]) < 1 then 1 else 0) >= 3;
