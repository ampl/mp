var x{1..3} integer >= 0 <= 3;

minimize o: 6*x[1] + 5*x[2] - 3*x[3];

s.t. c1: (min(x[3], x[3]) >= 2) ==> (x[2] <= 0) else (x[2] <= 0);
s.t. c2: (min(x[3], x[3]) >= 2) or (x[2] + x[3] > 3);
s.t. c3: (x[2] + x[3] > 3) ==> (min(x[2], x[3]) <= 3);
s.t. c4: not (min(x[2], x[3]) <= 3);
