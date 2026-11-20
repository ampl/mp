var x{1..4} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (1)*x[2] + (6)*x[3] + (3)*x[4];
s.t. c1: (x[3] < 2) ==> (x[1] <= -2) else (x[3] <= -2);
s.t. c2: abs(x[1] - x[1]) >= 3;
s.t. c3: (if x[1] >= 2 then x[2] else x[1]) < 1;
s.t. c4: (x[3] < 2) ==> (x[1] <= -2) else (x[3] <= -2);
s.t. c5: (if x[1] >= 2 then x[2] else x[1]) < 1;
