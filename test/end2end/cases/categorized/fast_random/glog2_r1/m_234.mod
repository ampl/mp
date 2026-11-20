var x{1..3} integer >= -2 <= 3;
minimize o: (4)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: ((if x[3] != 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[3] < 3) <==> (x[2] + x[2] >= -1);
s.t. c3: numberof -2 in (x[1], x[2], x[3]) != 0;
