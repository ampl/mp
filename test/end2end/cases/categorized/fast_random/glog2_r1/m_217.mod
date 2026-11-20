var x{1..3} integer >= -2 <= 2;

maximize o: 1000*x[1] - 1000*x[2] + 17*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] <= 1) or (x[2] != 1) or (x[3] > 2);
s.t. c3: not ((x[3] != -1) and (x[1] != -1));
