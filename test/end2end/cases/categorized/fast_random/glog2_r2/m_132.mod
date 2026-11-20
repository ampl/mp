var x{1..4} integer >= -2 <= 1;
minimize o: (4)*x[1] + (5)*x[2] + (-5)*x[3] + (-1)*x[4];
s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4]) > 2;
s.t. c2: ((if x[2] == -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c3: (x[2] != -1) <==> (x[2] + x[1] <= -1);
s.t. c4: ((if x[4] != 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
