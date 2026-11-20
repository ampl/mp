var x{1..5} integer >= 0 <= 2;

maximize o: 2000*x[1] + 1000*x[2] - 666*x[3] - 3000*x[4] + 17*x[5];

s.t. c1: (x[2] <= 0) and (x[5] < 1) ==> (x[5] >= 0);
s.t. c2: ((if x[5] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
