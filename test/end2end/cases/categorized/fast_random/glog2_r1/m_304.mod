var x{1..4} integer >= -2 <= 2;

maximize o: x[1] - x[2] - x[3] - 3*x[4];

s.t. c1: not ((x[4] != 2) and (x[4] < -1));
s.t. c2: exists {i in 1..4} x[i] == 2;
