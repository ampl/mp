var x{1..6} integer >= -2 <= 1;
minimize o: (-666)*x[1] + (34)*x[2] + (-666)*x[3] + (-2000)*x[4] + (-3000)*x[5] + (-333)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= -2);
s.t. c2: alldiff(x[5], x[1], x[4], x[6]);
s.t. c3: alldiff(x[5], x[6]);
