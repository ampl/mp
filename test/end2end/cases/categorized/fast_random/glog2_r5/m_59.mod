var x{1..4} integer >= 0 <= 4;
maximize o: (-2)*x[1] + (-5)*x[2] + (1)*x[3] + (0)*x[4];
s.t. c1: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c2: (x[1] <= 1) ==> (x[3] <= 4) else (x[4] <= 0);
s.t. c3: ((if x[2] >= 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c4: abs(x[3] - x[3]) != 4;
s.t. c5: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c6: abs(x[3] - x[3]) != 4;
