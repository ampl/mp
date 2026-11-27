var x{1..4} integer >= -2 <= 2;

maximize o: -x[1] + 2*x[2] - 2*x[3] - 4*x[4];

s.t. c1: ((if x[3] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[1] >= -2) and (x[1] >= 0) ==> (x[4] >= -2);
