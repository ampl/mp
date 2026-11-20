var x{1..6} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (4)*x[2] + (-2)*x[3] + (4)*x[4] + (2)*x[5] + (2)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= 1);
s.t. c2: (x[5] >= 1) or (x[5] < 1) or (x[5] > 3);
s.t. c3: abs(x[2] - x[1]) > 0;
s.t. c4: (x[6] > 1) ==> (x[1] >= 2) else (x[4] <= 0);
