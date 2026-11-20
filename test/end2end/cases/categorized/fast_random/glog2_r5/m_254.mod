var x{1..4} integer >= -3 <= 2;

minimize o: x[1] + 6*x[2] + x[3] - x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: abs(x[2] - x[1]) > 5;
s.t. c3: not ((x[3] >= -3) and (x[1] != 0));
s.t. c4: (x[2] != -3) and (x[2] != -2) ==> (x[2] >= -3);
