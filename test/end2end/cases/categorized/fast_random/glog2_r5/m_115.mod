var x{1..5} integer >= 0 <= 2;

minimize o: -3*x[1] + 6*x[2] - 5*x[3] + 3*x[4] + 0*x[5];

s.t. c1: (x[1] > 1) ==> (x[4] >= 0) else (x[2] >= 2);
s.t. c2: (x[5] > 2) and (x[4] == 2) ==> (x[1] >= 2);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) != 0;
s.t. c4: alldiff(x[1], x[3], x[5], x[2]);
s.t. c5: atmost 1 {i in 1..5} (x[i] >= 2);
s.t. c6: (x[4] > 2) or (x[2] < 2) or (x[4] > 1);
