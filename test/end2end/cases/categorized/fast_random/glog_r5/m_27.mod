var x{1..5} integer >= 0 <= 2;

minimize o: 5*x[1] - 4*x[2] + 0*x[3] + 6*x[4] - 5*x[5];

s.t. c1: (x[1] <= 0) or (x[2] <= 1) or (x[2] != 0);
s.t. c2: ((if x[1] <= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[4] == 1) <==> (x[2] + x[3] < 2);
s.t. c4: (if x[4] == 0 then x[3] else x[4]) <= 1;
