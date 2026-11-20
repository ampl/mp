var x{1..3} integer >= 0 <= 3;
maximize o: (1)*x[1] + (-5)*x[2] + (6)*x[3];
s.t. c1: abs(x[3] - x[2]) == 2;
s.t. c2: (x[1] <= 3) or (x[3] > 1) or (x[3] <= 2);
s.t. c3: ((if x[2] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
