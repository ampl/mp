var x{1..5} integer >= 0 <= 4;

minimize o: -5*x[1] - x[2] + 4*x[3] + 6*x[4] - 3*x[5];

s.t. c1: not ((x[1] <= 1) and (x[1] != 1));
s.t. c2: alldiff(x[2], x[3], x[4], x[1]);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c4: (x[4] != 0) <==> (x[2] + x[4] <= 4);
s.t. c5: exists {i in 1..5} x[i] == 2;
