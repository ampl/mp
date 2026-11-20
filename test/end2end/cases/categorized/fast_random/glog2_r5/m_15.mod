var x{1..6} integer >= -2 <= 0;

maximize o: 5*x[1] - 2*x[2] + 0*x[3] - 2*x[4] + 4*x[5] - 5*x[6];

s.t. c1: (x[1] > -2) ==> (x[3] == -2) else (x[3] >= 0);
s.t. c2: not ((x[2] < -1) and (x[5] > -1));
s.t. c3: max(x[5], x[5]) - min(x[6], x[6]) <= 0;
s.t. c4: alldiff(x[5], x[6]);
s.t. c5: numberof -2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
