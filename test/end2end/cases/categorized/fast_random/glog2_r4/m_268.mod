var x{1..7} integer >= 0 <= 1;

maximize o: -5*x[1] + 5*x[2] + 2*x[3] - 5*x[4] + x[5] + 6*x[6] + 4*x[7];

s.t. c1: (if x[3] >= 0 then x[6] else x[6]) < 0;
s.t. c2: ((if x[1] == 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[2] > 0) or (x[3] == 0) or (x[5] > 0);
s.t. c4: (x[1] >= 1) <==> (x[1] + x[4] == 1);
