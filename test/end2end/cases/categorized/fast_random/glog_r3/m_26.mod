var x{1..6} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (2)*x[2] + (-3)*x[3] + (-3)*x[4] + (-2)*x[5] + (6)*x[6];
s.t. c1: ((if x[3] <= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c2: (x[1] == 1) and (x[3] >= 0) ==> (x[1] >= 1);
s.t. c3: atleast 3 {i in 1..6} (x[i] >= -2);
