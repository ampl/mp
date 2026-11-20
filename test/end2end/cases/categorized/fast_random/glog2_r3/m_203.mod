var x{1..6} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (4)*x[2] + (4)*x[3] + (-3)*x[4] + (2)*x[5] + (-2)*x[6];
s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 0;
s.t. c2: ((if x[4] != 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: ((x[1] != -1) ==> (x[1] != -2)) and ((x[1] != -2) ==> (x[5] >= -1));
s.t. c4: exists {i in 1..6} x[i] == 2;
s.t. c5: x[6] == 2;
s.t. c6: x[2] == 1;
