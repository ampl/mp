var x{1..4} integer >= -2 <= 0;

minimize o: 0*x[1] - 3*x[2] + 0*x[3] + 2*x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= -2);
s.t. c2: (x[3] < 0) and (x[2] > 0) ==> (x[4] >= 0);
s.t. c3: atmost 1 {i in 1..4} (x[i] >= -2);
s.t. c4: (x[3] < 0) and (x[2] > 0) ==> (x[4] >= 0);
