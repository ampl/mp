var x{1..4} integer >= 0 <= 1;
minimize o: (0)*x[1] + (-2)*x[2] + (4)*x[3] + (5)*x[4];
s.t. c1: not ((x[4] == 1) and (x[3] <= 1));
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c3: (x[1] <= 1) <==> (x[3] + x[3] > 0);
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) >= 1;
