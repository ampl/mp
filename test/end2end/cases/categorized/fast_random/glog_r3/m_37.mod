var x{1..6} integer >= 0 <= 4;

maximize o: -x[1] + 4*x[2] - 4*x[3] - x[4] + x[5] + 0*x[6];

s.t. c1: ((if x[1] == 4 then 1 else 0) + (if x[4] >= 4 then 1 else 0)) >= 1;
s.t. c2: ((if x[6] != 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[1] <= 4) <==> (x[3] + x[2] <= 4);
