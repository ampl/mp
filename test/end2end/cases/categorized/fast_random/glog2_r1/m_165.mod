var x{1..4} integer >= -2 <= 1;
minimize o: (-5)*x[1] + (2)*x[2] + (-4)*x[3] + (-5)*x[4];
s.t. c1: abs(x[1] - x[4]) != 0;
s.t. c2: ((if x[2] != -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c3: max(x[3], x[4]) - min(x[1], x[1]) < -1;
s.t. c4: (x[2] != 1) or (x[4] != 0) or (x[1] < 0);
s.t. c5: atleast 1 {i in 1..4} (x[i] >= -1);
