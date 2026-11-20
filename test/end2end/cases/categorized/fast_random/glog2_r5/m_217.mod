var x{1..6} integer >= 0 <= 4;
minimize o: (51)*x[1] + (-2000)*x[2] + (-2000)*x[3] + (3000)*x[4] + (-333)*x[5] + (-666)*x[6];
s.t. c1: ((if x[6] < 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
