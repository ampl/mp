var x{1..4} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (4)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c2: (if x[2] != 0 then x[3] else x[1]) < 2;
s.t. c3: (x[2] >= 1) and (x[4] > 1) ==> (x[3] >= 1);
s.t. c4: x[1] == 2;
