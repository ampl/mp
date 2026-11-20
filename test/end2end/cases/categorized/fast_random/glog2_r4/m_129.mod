var x{1..5} integer >= -2 <= 0;

maximize o: -3*x[1] + 2*x[2] + 6*x[3] + 5*x[4] + 4*x[5];

s.t. c1: ((if x[5] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: abs(x[5] - x[1]) < 0;
s.t. c3: ((if x[5] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: abs(x[5] - x[1]) < 0;
