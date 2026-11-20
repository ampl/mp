var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 4*x[2] + 4*x[3] + x[4];

s.t. c1: ((x[1] <= 0) ==> (x[3] <= 0)) and ((x[3] <= 0) ==> (x[3] >= 0));
s.t. c2: ((if x[1] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c3: x[2] == 1;
