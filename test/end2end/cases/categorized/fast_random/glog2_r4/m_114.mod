var x{1..3} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (1)*x[2] + (5)*x[3];
s.t. c1: ((if x[2] <= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[1] >= 2) <==> (x[3] + x[3] < 1);
