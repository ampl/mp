var x{1..4} integer >= 0 <= 2;

maximize o: 2*x[1] + 2*x[2] - 2*x[3] + 2*x[4];

s.t. c1: not ((x[3] < 1) and (x[2] != 0));
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: x[2] == 2;
