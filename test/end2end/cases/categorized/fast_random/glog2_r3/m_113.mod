var x{1..6} integer >= -2 <= 1;
minimize o: (-3)*x[1] + (1)*x[2] + (-4)*x[3] + (2)*x[4] + (4)*x[5] + (6)*x[6];
s.t. c1: (x[4] <= -1) <==> (x[1] + x[6] < 1);
s.t. c2: ((if x[4] < -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 3;
s.t. c4: (if x[3] >= 0 then x[6] else x[2]) < 0;
s.t. c5: x[3] == -1;
s.t. c6: x[5] == -2;
