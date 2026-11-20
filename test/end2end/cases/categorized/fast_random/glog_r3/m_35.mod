var x{1..5} integer >= -2 <= 0;

maximize o: -3*x[1] + 2*x[2] - 5*x[3] + 3*x[4] - x[5];

s.t. c1: not ((x[3] == -1) and (x[4] == -2));
s.t. c2: (x[2] >= 0) <==> (x[2] + x[3] >= -2);
s.t. c3: ((if x[3] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
