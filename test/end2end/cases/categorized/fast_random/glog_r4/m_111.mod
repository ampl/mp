var x{1..5} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-3)*x[2] + (3)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) != 0;
s.t. c2: (x[1] > 1) ==> (x[2] <= 0) else (x[4] <= 1);
