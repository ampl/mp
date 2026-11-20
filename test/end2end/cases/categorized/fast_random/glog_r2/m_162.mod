var x{1..3} integer >= 0 <= 3;

minimize o: x[1] + 0*x[2] + 6*x[3];

s.t. c1: alldiff(x[2], x[1]);
s.t. c2: (x[2] <= 0) <==> (x[1] + x[3] >= 3);
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 1);
