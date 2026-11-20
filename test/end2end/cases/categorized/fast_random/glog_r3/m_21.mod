var x{1..3} integer >= -2 <= 2;

minimize o: 5*x[1] - 5*x[2] - 3*x[3];

s.t. c1: max(x[3], x[1]) - min(x[3], x[2]) == 2;
s.t. c2: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[2] >= 1) <==> (x[2] + x[2] > -2);
