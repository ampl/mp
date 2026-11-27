var x{1..5} integer >= -2 <= 2;

maximize o: -2*x[1] - 2*x[2] - x[3] + 5*x[4] + 2*x[5];

s.t. c1: (x[2] > 1) <==> (x[5] + x[1] >= 2);
s.t. c2: ((if x[4] > -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
s.t. c3: (x[5] != -2) <==> (x[4] + x[4] >= -1);
s.t. c4: abs(x[3] - x[3]) >= 2;
