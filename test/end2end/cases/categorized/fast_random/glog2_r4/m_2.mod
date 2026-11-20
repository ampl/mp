var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] + 0*x[2] + 2*x[3] + x[4] + 4*x[5];

s.t. c1: not ((x[1] < 0) and (x[4] == 0));
s.t. c2: (x[1] == 3) ==> (x[4] < 2) else (x[1] >= 2);
s.t. c3: ((if x[4] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
