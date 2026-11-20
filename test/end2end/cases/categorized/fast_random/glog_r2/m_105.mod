var x{1..4} integer >= -2 <= 2;
minimize o: (5)*x[1] + (-1)*x[2] + (6)*x[3] + (4)*x[4];
s.t. c1: ((if x[1] >= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] < 1) or (x[2] < 2) or (x[3] <= -1);
s.t. c3: ((if x[2] < 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
