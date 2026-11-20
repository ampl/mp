var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + 3*x[2] + x[3] + 6*x[4];

s.t. c1: (x[2] == 1) ==> (x[3] > 2);
s.t. c2: exists {i in 1..4} x[i] == 1;
