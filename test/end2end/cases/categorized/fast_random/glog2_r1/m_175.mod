var x{1..5} integer >= -2 <= 0;
maximize o: (0)*x[1] + (4)*x[2] + (3)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: atleast 2 {i in 1..5} (x[i] >= 0);
s.t. c2: not ((x[3] <= -1) and (x[2] < -1));
s.t. c3: (x[3] < -2) ==> (x[1] != 0) else (x[5] == -1);
s.t. c4: ((if x[5] >= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c5: ((if x[1] == -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
