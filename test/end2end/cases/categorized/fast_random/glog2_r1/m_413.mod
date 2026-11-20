var x{1..3} integer >= 0 <= 3;

minimize o: 2*x[1] + 0*x[2] - 2*x[3];

s.t. c1: not ((x[2] == 0) and (x[3] < 1));
s.t. c2: not ((x[3] <= 0) and (x[2] != 2));
s.t. c3: exists {i in 1..3} x[i] == 0;
s.t. c4: x[2] == 3;
