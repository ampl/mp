var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] - 5*x[2] + 3*x[3] + 2*x[4];

s.t. c1: (x[3] > 2) ==> (x[3] <= 4) else (x[3] <= 3);
s.t. c2: exists {i in 1..4} x[i] == 2;
s.t. c3: x[2] == 0;
