var x{1..6} integer >= 0 <= 1;

maximize o: -4*x[1] - 5*x[2] - 5*x[3] - 4*x[4] - 5*x[5] + 4*x[6];

s.t. c1: (x[4] >= 1) <==> (x[5] + x[3] < 0);
s.t. c2: (x[3] <= 1) or (x[6] != 0) or (x[4] != 1);
s.t. c3: ((if x[1] != 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
