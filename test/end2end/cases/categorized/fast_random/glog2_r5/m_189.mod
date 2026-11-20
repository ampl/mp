var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] - 2*x[2] - 3*x[3] + 5*x[4];

s.t. c1: (x[3] <= 1) ==> (x[2] >= 1) else (x[2] == 0);
s.t. c2: max(x[4], x[3]) - min(x[2], x[2]) > 1;
s.t. c3: (x[3] <= 1) ==> (x[2] >= 1) else (x[2] == 0);
s.t. c4: max(x[4], x[3]) - min(x[2], x[2]) > 1;
