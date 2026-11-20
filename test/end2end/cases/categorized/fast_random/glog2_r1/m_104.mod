var x{1..4} integer >= -3 <= 1;
minimize o: (-3)*x[1] + (-1)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: (x[2] <= -1) ==> (x[4] == -1) else (x[1] < 0);
s.t. c2: alldiff(x[2], x[3], x[4]);
