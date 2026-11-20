var x{1..3} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-4)*x[2] + (6)*x[3];
s.t. c1: ((x[1] <= 3) ==> (x[1] < 0)) and ((x[1] < 0) ==> (x[1] >= 3));
s.t. c2: abs(x[1] - x[3]) >= 1;
s.t. c3: exists {i in 1..3} x[i] == 0;
s.t. c4: ((if x[1] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
