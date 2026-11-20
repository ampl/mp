var x{1..4} integer >= -2 <= 2;

minimize o: 0*x[1] - x[2] + 3*x[3] - 3*x[4];

s.t. c1: (x[2] <= 2) ==> (x[4] != -2) else (x[2] != -1);
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: (x[1] == 2) ==> (x[3] > -2) else (x[1] == 0);
s.t. c4: atmost 3 {i in 1..4} (x[i] >= 2);
