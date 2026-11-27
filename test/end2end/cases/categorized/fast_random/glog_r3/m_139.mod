var x{1..5} integer >= -2 <= 2;

minimize o: 5*x[1] - x[2] - x[3] - 5*x[4] + 6*x[5];

s.t. c1: max(x[4], x[3]) - min(x[3], x[1]) > 1;
s.t. c2: (x[5] != 2) and (x[5] <= 1) ==> (x[4] >= 2);
s.t. c3: alldiff(x[3], x[4], x[2], x[5]);
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 0;
