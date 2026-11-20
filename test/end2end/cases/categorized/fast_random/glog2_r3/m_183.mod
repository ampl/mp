var x{1..6} integer >= 0 <= 4;
minimize o: (5)*x[1] + (-2)*x[2] + (-3)*x[3] + (4)*x[4] + (-3)*x[5] + (4)*x[6];
s.t. c1: ((if x[1] < 3 then 1 else 0) + (if x[6] >= 3 then 1 else 0)) >= 1;
s.t. c2: ((if x[2] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c3: atleast 1 {i in 1..6} (x[i] >= 1);
s.t. c4: x[3] == 4;
