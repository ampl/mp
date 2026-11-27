var x{1..4} integer >= 0 <= 3;

minimize o: 4*x[1] - 2*x[2] + x[3] + 5*x[4];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) >= 1;
s.t. c2: max(x[3], x[3]) - min(x[1], x[1]) >= 0;
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) == 1;
s.t. c5: (x[1] <= 1) or (x[1] > 3) or (x[3] <= 0);
