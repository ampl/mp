var x{1..4} integer >= -2 <= 1;
maximize o: (1000)*x[1] + (3000)*x[2] + (-333)*x[3] + (1000)*x[4];
s.t. c1: ((if x[4] <= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c2: abs(x[1] - x[4]) == 1;
s.t. c3: max(x[2], x[3]) - min(x[4], x[3]) != 0;
s.t. c4: (x[3] > -2) and (x[1] == 1) ==> (x[3] >= -2);
