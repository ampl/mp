var x{1..5} integer >= 0 <= 4;

maximize o: -999*x[1] + 1000*x[2] + 2000*x[3] - 999*x[4] + 17*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) <= 0;
s.t. c3: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) > 1;
s.t. c4: alldiff(x[2], x[3]);
