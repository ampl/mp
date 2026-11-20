var x{1..4} integer >= 0 <= 3;

maximize o: -2*x[1] - 4*x[2] + 3*x[3] + x[4];

s.t. c1: abs(x[2] - x[1]) != 2;
s.t. c2: ((if x[2] >= 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c3: (x[3] <= 2) <==> (x[1] + x[1] >= 0);
