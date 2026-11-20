var x{1..6} integer >= 0 <= 2;

minimize o: 6*x[1] - 5*x[2] - 2*x[3] - x[4] + 2*x[5] - 5*x[6];

s.t. c1: max(x[4], x[3]) - min(x[6], x[1]) <= 0;
s.t. c2: (x[5] < 0) ==> (x[3] == 0) else (x[5] <= 1);
