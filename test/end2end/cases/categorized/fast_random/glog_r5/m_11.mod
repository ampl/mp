var x{1..5} integer >= 0 <= 4;
minimize o: (-3)*x[1] + (-3)*x[2] + (5)*x[3] + (-1)*x[4] + (-1)*x[5];
s.t. c1: (x[3] != 2) <==> (x[5] + x[2] > 0);
s.t. c2: (if x[5] < 4 then x[1] else x[3]) == 4;
s.t. c3: not ((x[3] != 3) and (x[5] <= 4));
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) == 0;
