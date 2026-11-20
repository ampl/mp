var x{1..4} integer >= 0 <= 1;

maximize o: -5*x[1] - 4*x[2] - 2*x[3] - 2*x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= 0);
s.t. c2: ((if x[2] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[3] >= 0) and (x[4] <= 1) ==> (x[3] >= 0);
