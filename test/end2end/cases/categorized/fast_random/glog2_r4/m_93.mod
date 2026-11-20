var x{1..5} integer >= 0 <= 2;
minimize o: (1)*x[1] + (0)*x[2] + (0)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: alldiff(x[4], x[2]);
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 0);
s.t. c3: x[3] == 2;
s.t. c4: x[4] == 2;
