var x{1..4} integer >= -2 <= 0;
maximize o: (-3)*x[1] + (-2)*x[2] + (5)*x[3] + (1)*x[4];
s.t. c1: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c2: ((if x[4] > -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
s.t. c3: (x[4] > -2) ==> (x[4] < -2) else (x[4] <= -1);
s.t. c4: not ((x[4] <= 0) and (x[3] == -2));
