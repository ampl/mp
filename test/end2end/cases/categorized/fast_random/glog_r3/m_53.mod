var x{1..5} integer >= 0 <= 3;

maximize o: 4*x[1] + 4*x[2] - 5*x[3] + 3*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: ((if x[5] < 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((x[5] < 3) ==> (x[1] > 3)) and ((x[1] > 3) ==> (x[5] >= 3));
s.t. c4: alldiff(x[2], x[5], x[4]);
