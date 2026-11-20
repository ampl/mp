var x{1..5} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (-2)*x[2] + (1)*x[3] + (3)*x[4] + (-5)*x[5];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) <= 0;
s.t. c2: (if x[1] >= 0 then x[5] else x[4]) < 0;
s.t. c3: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c4: (if x[5] != 1 then x[2] else x[1]) != 2;
