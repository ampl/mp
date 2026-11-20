var x{1..5} integer >= 0 <= 4;

minimize o: -4*x[1] + 5*x[2] + 4*x[3] - 4*x[4] + 2*x[5];

s.t. c1: (x[1] > 1) ==> (x[1] <= 4);
s.t. c2: (x[3] == 3) and (x[4] != 4) ==> (x[4] >= 3);
s.t. c3: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) <= 1;
s.t. c4: alldiff(x[1], x[2], x[5], x[3]);
