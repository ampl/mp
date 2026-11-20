var x{1..3} integer >= -2 <= 0;

maximize o: 0*x[1] - x[2] - 4*x[3];

s.t. c1: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[1] == -2) <==> (x[1] + x[2] <= -2);
s.t. c3: exists {i in 1..3} x[i] == 0;
