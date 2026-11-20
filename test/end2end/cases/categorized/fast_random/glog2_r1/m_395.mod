var x{1..5} integer >= 0 <= 3;

maximize o: -2*x[1] - 3*x[2] + 6*x[3] + 2*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[5] <= 0) and (x[4] == 2) ==> (x[3] >= 0);
s.t. c3: not ((x[4] < 1) and (x[3] > 1));
s.t. c4: (x[4] == 2) or (x[4] >= 1) or (x[1] <= 0);
s.t. c5: (x[4] < 2) ==> (x[5] != 2) else (x[4] <= 1);
s.t. c6: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c7: (x[4] <= 1) ==> (x[1] < 1) else (x[4] < 0);
