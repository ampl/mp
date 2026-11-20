var x{1..4} integer >= 0 <= 3;

maximize o: -2*x[1] - 2*x[2] + 0*x[3] + 5*x[4];

s.t. c1: ((if x[2] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[4] >= 3) <==> (x[2] + x[3] > 1);
s.t. c3: x[4] == 1;
s.t. c4: x[3] == 2;
