var x{1..5} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (6)*x[2] + (2)*x[3] + (3)*x[4] + (-4)*x[5];
s.t. c1: ((if x[5] == 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 2;
s.t. c3: ((if x[5] == 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[2] != 1) <==> (x[5] + x[5] == 2);
s.t. c5: x[1] == 0;
