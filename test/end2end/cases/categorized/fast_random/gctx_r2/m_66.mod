var x{1..4} integer >= 0 <= 3;
minimize o: (1)*x[1] + (-1)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: (x[2] + x[4] > 1) ==> ((max(x[3], x[4]) != 2) or (x[1] <= 3));
s.t. c2: (if x[2] + x[4] > 1 then 1 else 0) + (if x[1] < 2 then 1 else 0) + (if x[1] < 2 then 1 else 0) >= 3;
s.t. c3: (x[1] < 2) or (x[2] + x[4] > 1);
s.t. c4: (if x[1] < 2 then 1 else 0) + (if x[2] + x[4] > 1 then 1 else 0) + (if max(x[3], x[4]) != 2 then 1 else 0) >= 2;
