var x{1..5} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: ((if x[1] <= -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
s.t. c2: ((x[2] < 2) ==> (x[1] == 1)) and ((x[1] == 1) ==> (x[3] >= 2));
