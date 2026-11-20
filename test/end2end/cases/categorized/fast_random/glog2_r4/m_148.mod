var x{1..6} integer >= 0 <= 1;

minimize o: 3*x[1] - 5*x[2] - 5*x[3] + 6*x[4] - 2*x[5] + x[6];

s.t. c1: not ((x[5] == 1) and (x[2] < 0));
s.t. c2: (x[6] > 1) or (x[4] > 0) or (x[3] != 1);
s.t. c3: max(x[4], x[1]) - min(x[6], x[4]) >= 0;
s.t. c4: (x[4] >= 0) ==> (x[1] == 0) else (x[3] > 1);
