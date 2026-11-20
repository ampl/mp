var x{1..6} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (-3)*x[2] + (0)*x[3] + (5)*x[4] + (-5)*x[5] + (0)*x[6];
s.t. c1: ((if x[3] <= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] <= 1) <==> (x[4] + x[5] >= 0);
s.t. c3: abs(x[4] - x[2]) >= 2;
