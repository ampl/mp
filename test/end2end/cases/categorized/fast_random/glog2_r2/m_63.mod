var x{1..3} integer >= -2 <= 1;
minimize o: (6)*x[1] + (-3)*x[2] + (-1)*x[3];
s.t. c1: (if x[2] != -1 then x[2] else x[2]) != 0;
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: ((if x[1] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[3] == -2) ==> (x[1] <= -1) else (x[3] > -2);
s.t. c5: x[3] == -2;
