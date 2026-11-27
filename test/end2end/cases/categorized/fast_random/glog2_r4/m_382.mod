var x{1..3} integer >= 0 <= 4;

minimize o: -2*x[1] + 0*x[2] + 6*x[3];

s.t. c1: atmost 3 {i in 1..3} (x[i] >= 2);
s.t. c2: (x[3] <= 4) <==> (x[3] + x[3] >= 2);
s.t. c3: (x[3] <= 4) <==> (x[2] + x[2] > 0);
