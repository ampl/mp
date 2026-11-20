var x{1..6} integer >= 0 <= 4;
maximize o: (-666)*x[1] + (-3000)*x[2] + (-2000)*x[3] + (3000)*x[4] + (3000)*x[5] + (17)*x[6];
s.t. c1: max(x[3], x[4]) - min(x[5], x[5]) < 0;
s.t. c2: ((if x[4] > 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
