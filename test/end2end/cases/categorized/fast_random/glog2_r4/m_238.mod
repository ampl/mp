var x{1..5} integer >= 0 <= 1;

maximize o: 5*x[1] + 4*x[2] - x[3] - 3*x[4] + 6*x[5];

s.t. c1: ((if x[2] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: (x[3] <= 0) ==> (x[3] != 1);
s.t. c3: max(x[5], x[5]) - min(x[3], x[5]) != 0;
