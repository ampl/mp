var x{1..3} integer >= 0 <= 4;

maximize o: -3*x[1] + 2*x[2] - 3*x[3];

s.t. c1: (x[1] != 1) ==> (x[3] >= 4) else (x[3] < 2);
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (x[2] > 3) ==> (x[1] > 0);
