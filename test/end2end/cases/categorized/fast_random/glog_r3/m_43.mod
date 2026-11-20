var x{1..6} integer >= -2 <= 1;
maximize o: (-3)*x[1] + (4)*x[2] + (0)*x[3] + (-2)*x[4] + (1)*x[5] + (6)*x[6];
s.t. c1: (x[3] <= -1) ==> (x[6] == -1);
s.t. c2: (x[1] < -2) or (x[6] != 1) or (x[5] > -2);
s.t. c3: ((if x[6] <= -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
s.t. c4: (if x[2] < 0 then x[1] else x[6]) < -1;
