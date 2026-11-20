var x{1..6} integer >= 0 <= 2;
maximize o: (34)*x[1] + (51)*x[2] + (-999)*x[3] + (-666)*x[4] + (-333)*x[5] + (1000)*x[6];
s.t. c1: ((if x[1] != 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: atleast 2 {i in 1..6} (x[i] >= 2);
