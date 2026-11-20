var x{1..6} integer >= -2 <= 1;

maximize o: 2*x[1] - 5*x[2] + 6*x[3] - 2*x[4] + 4*x[5] + 4*x[6];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
s.t. c2: not ((x[4] >= 0) and (x[1] > -1));
s.t. c3: (x[1] <= -1) and (x[4] >= 1) ==> (x[4] >= -1);
