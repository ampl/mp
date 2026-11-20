var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[2] < 1) and (x[2] >= 1) ==> (x[4] >= 1);
s.t. c2: (if x[4] >= 1 then x[1] else x[3]) == 0;
s.t. c3: ((if x[2] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
