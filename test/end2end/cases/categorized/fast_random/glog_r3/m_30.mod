var x{1..5} integer >= -2 <= 0;
minimize o: (2)*x[1] + (-3)*x[2] + (-1)*x[3] + (1)*x[4] + (4)*x[5];
s.t. c1: not ((x[4] == 0) and (x[1] <= -1));
s.t. c2: (x[1] != 0) or (x[5] <= -2) or (x[1] == 0);
s.t. c3: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c4: (x[4] != 0) ==> (x[5] >= -1) else (x[5] != -2);
