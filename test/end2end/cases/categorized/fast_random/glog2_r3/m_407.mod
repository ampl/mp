var x{1..6} integer >= 0 <= 2;
minimize o: (2000)*x[1] + (51)*x[2] + (-666)*x[3] + (2000)*x[4] + (3000)*x[5] + (17)*x[6];
s.t. c1: ((if x[3] < 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..6} x[i] == 1;
