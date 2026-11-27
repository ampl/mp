var x{1..4} integer >= -2 <= 2;

minimize o: 6*x[1] + 4*x[2] + 4*x[3] - 5*x[4];

s.t. c1: atmost 2 {i in 1..4} (x[i] >= 2);
s.t. c2: ((if x[4] >= -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[3], x[4]);
s.t. c4: (x[1] != 2) and (x[3] > 0) ==> (x[2] >= 2);
