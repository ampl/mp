var x{1..3} integer >= -2 <= 2;

minimize o: 4*x[1] - 2*x[2] + 3*x[3];

s.t. c1: atmost 1 {i in 1..3} (x[i] >= -2);
s.t. c2: (x[2] > 3) ==> (x[2] != 3) else (x[1] > -2);
s.t. c3: (x[3] > -2) ==> (x[3] > 3) else (x[2] <= 2);
