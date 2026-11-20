var x{1..4} integer >= 0 <= 3;
minimize o: (0)*x[1] + (2)*x[2] + (1)*x[3] + (4)*x[4];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) != 0;
s.t. c2: (x[3] <= 3) ==> (x[2] == 3) else (x[1] <= 0);
s.t. c3: ((x[3] > 2) ==> (x[3] > 0)) and ((x[3] > 0) ==> (x[4] >= 2));
s.t. c4: alldiff(x[1], x[2], x[4]);
