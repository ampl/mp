var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] + 6*x[2] + 6*x[3];

s.t. c1: (x[3] >= 1) ==> (x[1] > 4) else (x[3] >= 0);
s.t. c2: not ((x[2] != 0) and (x[2] >= 3));
s.t. c3: abs(x[2] - x[3]) > 4;
