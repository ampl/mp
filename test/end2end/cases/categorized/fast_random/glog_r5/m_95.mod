var x{1..6} integer >= 0 <= 2;

maximize o: -2*x[1] + 6*x[2] + 2*x[3] - 2*x[4] - 5*x[5] + 5*x[6];

s.t. c1: max(x[2], x[5]) - min(x[4], x[2]) >= 2;
s.t. c2: (x[1] <= 0) ==> (x[4] != 2) else (x[6] == 2);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
