var x{1..3} integer >= 0 <= 4;

maximize o: 4*x[1] - 4*x[2] + 5*x[3];

s.t. c1: (x[1] != 1) <==> (x[3] + x[2] < 3);
s.t. c2: ((if x[1] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: (if x[2] > 4 then x[2] else x[3]) <= 3;
