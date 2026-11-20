var x{1..4} integer >= 0 <= 4;

minimize o: 2*x[1] + 6*x[2] + 4*x[3] - 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == 3;
s.t. c2: not ((x[3] != 4) and (x[3] >= 2));
s.t. c3: (x[3] < 3) or (x[1] == 0) or (x[4] >= 4);
