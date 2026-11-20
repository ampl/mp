var x{1..4} integer >= -2 <= 2;

minimize o: 5*x[1] + x[2] - 2*x[3] + 6*x[4];

s.t. c1: (x[1] == 0) and (x[1] < 2) ==> (x[2] >= 0);
s.t. c2: not ((x[2] >= 2) and (x[3] >= -1));
s.t. c3: not ((x[2] > 0) and (x[3] == -1));
