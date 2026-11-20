var x{1..6} integer >= 0 <= 1;

maximize o: -4*x[1] - 2*x[2] + 2*x[3] - 2*x[4] - x[5] - 2*x[6];

s.t. c1: ((if x[4] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: ((x[6] >= 0) ==> (x[1] >= 0)) and ((x[1] >= 0) ==> (x[1] >= 0));
s.t. c3: (x[2] != 1) ==> (x[6] == 0) else (x[6] == 1);
