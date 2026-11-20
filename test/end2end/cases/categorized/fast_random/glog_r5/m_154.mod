var x{1..4} integer >= -2 <= 0;
minimize o: (-3)*x[1] + (1)*x[2] + (-1)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[2], x[3], x[4]);
s.t. c2: ((if x[2] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[1] != -2) or (x[4] == -1) or (x[3] < -1);
