var x{1..3} integer >= -1 <= 2;
maximize o: (2)*x[1] + (-5)*x[2] + (6)*x[3];
s.t. c1: (if max(x[3], x[2]) > -1 then 1 else 0) + (if x[3] + x[2] != 2 then 1 else 0) + (if x[3] + x[2] != 2 then 1 else 0) >= 3;
s.t. c2: (if max(x[3], x[2]) > -1 then 1 else 0) + (if x[3] + x[2] != 2 then 1 else 0) + (if x[1] + x[3] <= 2 then 1 else 0) >= 3;
s.t. c3: (x[3] + x[2] != 2) or (x[1] + x[3] <= 2);
