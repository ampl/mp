var x{1..6} integer >= -2 <= 1;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: atleast 2 {i in 1..6} (x[i] >= -1);
s.t. c2: alldiff(x[4], x[6]);
s.t. c3: (if x[4] >= 1 then x[5] else x[2]) > 1;
