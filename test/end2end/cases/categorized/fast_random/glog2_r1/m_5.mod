var x{1..6} integer >= 0 <= 3;
maximize o: (2)*x[1] + (6)*x[2] + (-2)*x[3] + (0)*x[4] + (-5)*x[5] + (2)*x[6];
s.t. c1: abs(x[2] - x[3]) < 0;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: ((if x[6] >= 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
s.t. c4: atleast 2 {i in 1..6} (x[i] >= 1);
s.t. c5: ((if x[2] <= 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
