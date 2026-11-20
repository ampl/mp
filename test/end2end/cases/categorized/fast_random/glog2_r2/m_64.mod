var x{1..5} integer >= -3 <= 1;

maximize o: 3*x[1] + 6*x[2] + 3*x[3] + 3*x[4] + 5*x[5];

s.t. c1: not ((x[2] > 1) and (x[3] == -1));
s.t. c2: (x[4] < -1) ==> (x[4] != -3);
s.t. c3: max(x[3], x[2]) - min(x[5], x[3]) > -1;
s.t. c4: ((if x[1] >= 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
