var x{1..5} integer >= -2 <= 1;
minimize o: (51)*x[1] + (-666)*x[2] + (17)*x[3] + (-333)*x[4] + (-1000)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[5], x[1]);
s.t. c3: ((x[2] <= -1) ==> (x[5] > 1)) and ((x[5] > 1) ==> (x[4] >= -1));
s.t. c4: ((if x[2] <= -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
