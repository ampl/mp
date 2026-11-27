var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (x[1] <= 1) <==> (x[3] + x[1] <= 0);
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: (if x[1] > 1 then x[2] else x[3]) <= 1;
