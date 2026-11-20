var x{1..4} integer >= 0 <= 4;

minimize o: 5*x[1] - x[2] + 0*x[3] + x[4];

s.t. c1: atmost 3 {i in 1..4} (x[i] >= 1);
s.t. c2: (x[2] < 4) and (x[1] >= 1) ==> (x[2] >= 4);
s.t. c3: not ((x[2] >= 4) and (x[3] <= 2));
s.t. c4: (x[4] > 0) and (x[3] == 2) ==> (x[1] >= 0);
s.t. c5: atmost 3 {i in 1..4} (x[i] >= 1);
s.t. c6: (x[4] > 0) and (x[3] == 2) ==> (x[1] >= 0);
