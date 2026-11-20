var x{1..5} integer >= -2 <= 0;

minimize o: -999*x[1] + 3000*x[2] + 1000*x[3] - 2000*x[4] - 3000*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[3], x[1]);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
