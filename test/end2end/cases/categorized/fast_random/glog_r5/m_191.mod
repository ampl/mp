var x{1..3} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (-5)*x[2] + (-5)*x[3];
s.t. c1: not ((x[1] == 1) and (x[3] > 0));
s.t. c2: numberof -1 in (x[1], x[2], x[3]) != 1;
s.t. c3: alldiff(x[2], x[1], x[3]);
s.t. c4: (x[3] >= -1) ==> (x[2] < -2) else (x[2] == -1);
