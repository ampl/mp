var x{1..6} integer >= 0 <= 2;

minimize o: 2*x[1] + 5*x[2] + 3*x[3] + 5*x[4] + 6*x[5] + 6*x[6];

s.t. c1: max(x[5], x[5]) - min(x[5], x[5]) != 2;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 0;
s.t. c3: (x[5] > 0) ==> (x[5] == 2) else (x[5] > 0);
s.t. c4: abs(x[5] - x[3]) <= 2;
