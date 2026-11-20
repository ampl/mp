var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] - x[2] + 4*x[3] - 3*x[4];

s.t. c1: (x[3] <= 0) and (x[4] != 1) ==> (x[3] >= 0);
s.t. c2: not ((x[2] >= 4) and (x[3] < 0));
s.t. c3: exists {i in 1..4} x[i] == 0;
