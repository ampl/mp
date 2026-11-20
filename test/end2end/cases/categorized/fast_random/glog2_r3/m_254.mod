var x{1..5} integer >= -3 <= 1;
minimize o: (2)*x[1] + (4)*x[2] + (-3)*x[3] + (3)*x[4] + (1)*x[5];
s.t. c1: atleast 2 {i in 1..5} (x[i] >= 0);
s.t. c2: abs(x[5] - x[5]) > 0;
s.t. c3: (x[2] > 0) ==> (x[5] < 1) else (x[3] < -3);
