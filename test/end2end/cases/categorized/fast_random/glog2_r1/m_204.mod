var x{1..5} integer >= -2 <= 3;
maximize o: (3)*x[1] + (4)*x[2] + (1)*x[3] + (3)*x[4] + (0)*x[5];
s.t. c1: ((if x[3] <= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: not ((x[5] < -2) and (x[2] <= 3));
s.t. c3: (x[3] <= -1) ==> (x[5] < 3) else (x[2] > -1);
