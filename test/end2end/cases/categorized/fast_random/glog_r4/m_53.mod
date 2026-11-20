var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (3)*x[2] + (3)*x[3];
s.t. c1: (x[1] <= 2) ==> (x[3] > 1) else (x[2] > 2);
s.t. c2: max(x[1], x[3]) - min(x[1], x[1]) <= 1;
s.t. c3: not ((x[3] == 2) and (x[2] >= 2));
