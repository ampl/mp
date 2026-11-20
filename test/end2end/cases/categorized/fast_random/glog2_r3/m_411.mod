var x{1..5} integer >= 0 <= 4;

maximize o: -x[1] - 4*x[2] + 0*x[3] + 5*x[4] + 2*x[5];

s.t. c1: (x[5] == 0) ==> (x[5] == 4) else (x[2] == -1);
s.t. c2: ((if x[1] == 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
