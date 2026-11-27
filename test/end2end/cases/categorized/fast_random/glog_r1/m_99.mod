var x{1..3} integer >= -2 <= 0;

maximize o: 4*x[1] - x[2] + 3*x[3];

s.t. c1: (x[2] == 0) ==> (x[1] >= -1) else (x[3] > 0);
s.t. c2: exists {i in 1..3} x[i] == -1;
