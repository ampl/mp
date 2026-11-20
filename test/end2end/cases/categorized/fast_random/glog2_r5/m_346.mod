var x{1..6} integer >= -2 <= 0;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: abs(x[4] - x[6]) >= 1;
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: ((if x[4] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c4: (x[2] <= -2) ==> (x[1] < -2);
