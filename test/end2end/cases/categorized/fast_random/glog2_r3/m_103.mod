var x{1..6} integer >= -2 <= 0;
maximize o: (-2)*x[1] + (-4)*x[2] + (3)*x[3] + (6)*x[4] + (2)*x[5] + (5)*x[6];
s.t. c1: ((if x[2] < 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
s.t. c2: abs(x[4] - x[2]) == 2;
s.t. c3: atleast 2 {i in 1..6} (x[i] >= 0);
s.t. c4: (x[6] != 0) and (x[3] <= -2) ==> (x[4] >= 0);
s.t. c5: x[1] == -2;
